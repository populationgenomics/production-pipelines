import re
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.jobs.samtools import samtools_stats

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='samtools-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'samtools': 'test_image',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
            },
        },
    )


def setup_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    batch = create_local_batch(tmp_path)
    cram_pth = create_cram_input(location=tmp_path, index=True)

    return config, cram_pth, batch


class TestSamtoolsStatsRun:
    def test_creates_one_job(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        _ = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
            overwrite=False,
        )

        assert (
            len(batch.select_jobs('samtools stats')) == 1
        ), 'Unexpected number of samtools stats jobs in batch list, should be just 1 job'

    def test_will_return_none_if_path_already_exists(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_stats_file'
        output.touch()

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=False,
        )

        assert j is None

    def test_will_create_job_if_path_already_exists_and_overwrite_true(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)
        output = tmp_path / 'output_stats_file'
        output.touch()

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=True,
        )

        assert j is not None

    @pytest.mark.parametrize(
        'job_attrs, expected_attrs',
        [
            (None, {'tool': 'samtools'}),
            (
                {'test_tool': 'test_samtools'},
                {'test_tool': 'test_samtools', 'tool': 'samtools'},
            ),
        ],
    )
    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path, job_attrs, expected_attrs):
        _, cram_pth, batch = setup_test(tmp_path)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=job_attrs,
            overwrite=False,
        )

        assert j is not None
        assert j.attributes == expected_attrs

    def test_uses_samtools_image_specified_in_config(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        assert j is not None
        assert j._image == config.images['samtools']

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        config = default_config()
        config.workflow.ref_fasta = 'test_workflow_ref.fa'
        _, cram_pth, batch = setup_test(tmp_path, config)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'--reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        ref_file = config.references['broad']['ref_fasta']
        assert re.search(fr'--reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        # This test ensures the script uses a fail-safe copy operation when working with CRAM files and their indices.
        # Such retries are common for network calls like copying files from GCS to handle unexpected issues.
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    def test_samtools_writes_to_resource_file_named_output_stats(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        assert re.search(r'\$CRAM > \${BATCH_TMPDIR}/.*/output_stats', cmd)

    def test_batch_writes_samtools_stats_file_to_output_path(self, mocker: MockFixture, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=out_path,
            job_attrs=None,
        )

        assert j is not None
        spy.assert_called_with(j.output_stats, str(out_path))
