import re
from typing import Literal

import pandas as pd
import pytest
from cpg_utils import Path
from hailtop.batch import Batch, Resource
from hailtop.batch.job import Job
from pytest_mock import MockFixture

from cpg_workflows.jobs.somalier import _check_pedigree, _relate, extract

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='somalier-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'somalier': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                },
                'somalier_sites': 'test_somalier_sites',
            },
            'cramqc': {'num_pcs': '4'},
        },
    )


def setup_test(
    tmp_path: Path,
    ref_fasta: str | None = None,
    sequencing_type: Literal['genome', 'exome'] | None = None,
    index: bool = True,
):
    config = default_config()

    if ref_fasta is not None:
        config.workflow.ref_fasta = ref_fasta

    if sequencing_type is not None:
        config.workflow.sequencing_type = sequencing_type

    set_config(config, tmp_path / 'config.toml')

    cram_pth = create_cram_input(
        location=tmp_path, prefix='test', index=index, reference_assembly='GRCh38.fa'
    )

    batch = create_local_batch(tmp_path)

    return config, cram_pth, batch


def setup_relate_test(tmp_path: Path):
    config = default_config()
    set_config(config, tmp_path / 'config.toml')

    dataset_id = config.workflow.dataset
    batch = create_local_batch(tmp_path)
    sg = create_sequencing_group(
        dataset=dataset_id,
        sequencing_type=config.workflow.sequencing_type,
    )
    somalier_path_by_sgid = {sg.id: (tmp_path / 'test.somalier')}

    return config, batch, sg, somalier_path_by_sgid


class TestSomalierExtract:
    def test_creates_one_job(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        _ = extract(b=batch, cram_path=cram_pth, out_somalier_path=tmp_path)

        # ---- Assertions
        assert (
            len(batch.select_jobs('Somalier extract')) == 1
        ), "Unexpected number of 'Somalier' jobs in batch list, should be just 1 job"

    def test_will_return_none_if_output_file_exists_and_overwrite_is_false(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=False,
        )

        # --- Assertions
        assert j is None, 'Job was created when output file should have been reused'

    def test_will_create_one_job_if_output_file_exists_and_overwrite_is_true(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=True,
        )

        # ---- Assertions
        assert len(batch.select_jobs('Somalier extract')) == 1
        assert j is not None

    def test_sets_supplied_job_attrs_or_sets_default_attrs_if_attrs_not_supplied(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j_default_attrs = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )
        j_supplied_attrs = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'second_output_file'),
            job_attrs={'test_tool': 'test_somalier'},
        )

        # ---- Assertions
        assert j_default_attrs is not None and j_supplied_attrs is not None
        assert j_default_attrs.attributes == {'tool': 'somalier'}
        assert j_supplied_attrs.attributes == {
            'test_tool': 'test_somalier',
            'tool': 'somalier',
        }

    def test_uses_image_specified_in_config(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        assert j is not None
        assert j._image == config.images['somalier']

    def test_sets_sites_location_with_name_of_sites_file_in_config(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        sites = config.other['references']['somalier_sites']
        assert j is not None
        assert re.search(fr'SITES=\$BATCH_TMPDIR/sites/{sites}', cmd)

    def test_runs_with_somalier_sites_from_config_file(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        sites = config.other['references']['somalier_sites']
        assert j is not None
        assert re.search(fr'--sites \${{BATCH_TMPDIR}}/inputs/\w+/{sites}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path, ref_fasta='test_workflow_ref.fa')

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        ref_file = config.other['references']['broad']['ref_fasta']
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_batch_writes_output_file_to_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=out_path,
        )

        # ---- Assertions
        assert j is not None
        spy.assert_called_with(j.output_file, str(out_path))

    def test_raises_error_if_no_cram_index_path_given(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path, index=False)

        # ---- The jobs we want to test and Assertions
        with pytest.raises(
            ValueError, match='CRAM for somalier is required to have CRAI index'
        ):
            _ = extract(
                b=batch,
                cram_path=cram_pth,
                out_somalier_path=(tmp_path / 'output_file'),
            )


class TestSomalierRelate:
    def test_creates_one_relate_job(self, tmp_path: Path):
        # ---- Test setup
        _, batch, sg, somalier_path_by_sgid = setup_relate_test(tmp_path)

        # ---- The job that we want to test
        j = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=None,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )

        # ---- Assertions
        relate_jobs = batch.select_jobs(rf'Somalier relate')
        assert len(relate_jobs) == 1

    def test_adds_label_to_default_job_title_if_label_provided(self, tmp_path: Path):
        # ---- Test setup
        _, batch, sg, somalier_path_by_sgid = setup_relate_test(tmp_path)

        label = 'test-label'

        # ---- The job that we want to test
        j = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=label,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )

        # ---- Assertions
        job_name = f'Somalier relate [{label}]'
        assert j.name == job_name

    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path):
        # ---- Test setup
        _, batch, sg, somalier_path_by_sgid = setup_relate_test(tmp_path)

        # ---- The job that we want to test
        j_default_attrs = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=None,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            job_attrs=None,
        )

        j_supplied_attrs = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=None,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            job_attrs={'test_tool': 'test_relate'},
        )

        # ---- Assertions
        assert j_default_attrs is not None and j_supplied_attrs is not None
        assert j_default_attrs.attributes == {'tool': 'somalier'}
        assert j_supplied_attrs.attributes == {
            'test_tool': 'test_relate',
            'tool': 'somalier',
        }

    def test_uses_image_specified_in_config(self, tmp_path: Path):
        # ---- Test setup
        config, batch, sg, somalier_path_by_sgid = setup_relate_test(tmp_path)

        # ---- The job that we want to test
        j = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=None,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )

        # ---- Assertions
        assert j is not None
        assert j._image == config.images['somalier']

    def test_if_verifybamid_exists_for_sg_check_freemix(self, tmp_path: Path):
        # ---- Test setup
        config, batch, sg, somalier_path_by_sgid = setup_relate_test(tmp_path)

        # ---- The job that we want to test
        j = _relate(
            b=batch,
            somalier_path_by_sgid=somalier_path_by_sgid,
            sequencing_group_ids=[sg.id],
            rich_id_map={sg.id: sg.pedigree.fam_id},
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            label=None,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            verifybamid_by_sgid=somalier_path_by_sgid,
        )

        # ---- Assertions
        cmd = get_command_str(j)
        verifybamid_file = somalier_path_by_sgid
        print()
        assert re.search(fr'FREEMIX=$(cat ${{BATCH_TMPDIR}}"/inputs/\w+')
