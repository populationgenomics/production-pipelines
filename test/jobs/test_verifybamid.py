import re
from typing import Literal

import pytest
from pytest_mock import MockFixture

from cpg_utils import Path
from cpg_workflows.jobs.verifybamid import verifybamid

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='verifybamid-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'verifybamid': 'test_image',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                'genome_contam_ud': 'test_genome_ud.ud',
                'genome_contam_bed': 'test_genome_bed.bed',
                'genome_contam_mu': 'test__genome_mu.mu',
                'exome_contam_ud': 'test_exome_ud.ud',
                'exome_contam_bed': 'test_exome_bed.bed',
                'exome_contam_mu': 'test_exome_mu.mu',
            },
        },
        other={'cramqc': {'num_pcs': '4'}},
    )


def setup_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    batch = create_local_batch(tmp_path)
    cram_pth = create_cram_input(location=tmp_path, index=True)

    return config, cram_pth, batch


class TestVerifyBAMID:
    def test_creates_one_job(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        _ = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        assert (
            len(batch.select_jobs('VerifyBamID')) == 1
        ), "Unexpected number of 'VerifyBamID' jobs in batch list, should be just 1 job"

    def test_will_return_none_if_path_already_exists(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_stats_file'
        output.touch()

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=output,
            job_attrs=None,
        )

        assert j is None, 'A new job was created when it should have reused existing output, was overwrite set to True?'

    def test_will_create_job_if_path_already_exists_and_overwrite_true(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)
        output = tmp_path / 'output_verifybamid_file'
        output.touch()

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=output,
            overwrite=True,
        )

        assert j is not None, 'Output not overwritten, no new job was created.'

    @pytest.mark.parametrize(
        'job_attrs, expected_attrs',
        [
            (None, {'tool': 'VerifyBamID'}),
            (
                {'test_tool': 'test_VerifyBamID'},
                {'test_tool': 'test_VerifyBamID', 'tool': 'VerifyBamID'},
            ),
        ],
    )
    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path, job_attrs, expected_attrs):
        _, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=job_attrs,
        )

        assert j is not None
        assert j.attributes == expected_attrs

    def test_uses_verifybamid_image_specified_in_config(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        assert j is not None
        assert j._image == config.images['verifybamid']

    def test_uses_num_principal_components_in_config_file_in_bash_command(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        num_pcs = config.other['cramqc']['num_pcs']
        cmd = get_command_str(j)
        assert re.search(fr'--NumPC {num_pcs}', cmd)

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        config = default_config()
        config.workflow.ref_fasta = 'test_workflow_ref.fa'
        _, cram_pth, batch = setup_test(tmp_path, config)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'--Reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        ref_file = config.references['broad']['ref_fasta']
        assert re.search(fr'--Reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_input_ud_bed_mu_files_set_in_config(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        sequencing_type = config.workflow.sequencing_type
        ud_file = config.references['broad'][f'{sequencing_type}_contam_ud']
        bed_file = config.references['broad'][f'{sequencing_type}_contam_bed']
        mu_file = config.references['broad'][f'{sequencing_type}_contam_mu']

        assert re.search(fr'--UDPath \${{BATCH_TMPDIR}}/inputs/\w+/{ud_file}', cmd)
        assert re.search(fr'--MeanPath \${{BATCH_TMPDIR}}/inputs/\w+/{mu_file}', cmd)
        assert re.search(fr'--BedPath \${{BATCH_TMPDIR}}/inputs/\w+/{bed_file}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        # This test ensures the script uses a fail-safe copy operation when working with CRAM files and their indices.
        # Such retries are common for network calls like copying files from GCS to handle unexpected issues.
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    @pytest.mark.parametrize('sequencing_type', ['exome', 'genome'])
    def test_extra_opts_changes_according_to_sequencing_type(
        self,
        tmp_path: Path,
        sequencing_type: Literal['exome', 'genome'],
    ):
        config = default_config()
        config.workflow.sequencing_type = sequencing_type
        config, cram_pth, batch = setup_test(tmp_path, config)

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        if sequencing_type == 'exome':
            assert re.search('--max-depth 1000', cmd)
        else:
            assert '--max-depth 1000' not in cmd

    def test_batch_writes_verifybamid_selfsm_file_to_output_path(self, mocker: MockFixture, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'

        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=out_path,
            job_attrs=None,
        )

        assert j is not None
        spy.assert_called_once_with(j.out_selfsm, str(out_path))
