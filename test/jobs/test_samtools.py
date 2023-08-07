import os
import re

from cpg_utils import Path
from cpg_utils.hail_batch import image_path
from pytest_mock import MockFixture

from cpg_workflows import utils
from cpg_workflows.batch import Batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs.samtools import samtools_stats

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
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
            }
        },
    )


def setup_test(tmp_path: Path, ref_fasta: str | None = None):
    config = default_config()

    if ref_fasta != None:
        config.workflow.ref_fasta = ref_fasta

    set_config(config, tmp_path / 'config.toml')

    cram_pth = create_cram_input(
        location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
    )

    batch = create_local_batch(tmp_path)

    return config, cram_pth, batch


class TestSamtoolsRun:
    def test_creates_one_job(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
            overwrite=False,
        )
        # ---- Assertions
        assert (
            len(batch.select_jobs('samtools stats')) == 1
        ), 'Unexpected number of samtools stats jobs in batch list, should be just 1 job'

    def test_will_return_none_if_path_already_exists(self, tmp_path: Path):
        # Giving output file that DOES exist
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_stats_file'
        output.touch()

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=False,
        )

        # --- Assertions
        assert j is None

    def test_will_create_job_if_path_already_exists_and_overwrite_true(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)
        output = tmp_path / 'output_stats_file'
        output.touch()

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=True,
        )

        # --- Assertions
        assert j is not None

    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j_default_attrs = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
            overwrite=False,
        )
        j_supplied_attrs = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'second_output_file'),
            job_attrs={'test_tool': 'test_samtools'},
            overwrite=False,
        )

        assert j_default_attrs.attributes == {'tool': 'samtools'}
        assert j_supplied_attrs.attributes == {
            'test_tool': 'test_samtools',
            'tool': 'samtools',
        }

    def test_uses_samtools_image_specified_in_config(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )
        assert j._image == config.images['samtools']

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        # Not sure how to check what reference file is being used. \
        # The resource group inside the job 'j' is not obvious in its naming

        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path, ref_fasta='test_workflow_ref.fa')

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'--reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        ref_file = config.references['broad']['ref_fasta']
        assert re.search(fr'--reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )

        cmd = get_command_str(j)
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    def test_samtools_writes_to_resource_file_named_output_stats(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_file'),
            job_attrs=None,
        )
        # ---- Assertions
        cmd = get_command_str(j)
        assert re.search(r'\$CRAM > \${BATCH_TMPDIR}/.*/output_stats', cmd)

    # TODO: What's the difference between this test and testing whether output exists?
    def test_batch_writes_samtools_stats_file_to_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=out_path,
            job_attrs=None,
        )

        spy.assert_called_with(j.output_stats, str(out_path))
