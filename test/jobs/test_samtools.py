import re
import os

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str

from cpg_workflows.jobs.samtools import samtools_stats
from cpg_utils import Path
from cpg_workflows.filetypes import CramPath
from cpg_utils.hail_batch import image_path


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
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                }
            },
        },
    )


class TestSamtoolsRun:
    def test_create_cram_input_method_creates_inputs(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )
        # ---- Check CRAM and CRAI files exist
        assert os.path.exists(tmp_path / 'test.cram')
        assert os.path.exists(tmp_path / 'test.cram.crai')

    def test_run_samtools(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )

        # ---- Assertions
        assert (
            len(batch.select_jobs('samtools stats')) == 1
        ), "Unexpected number of 'samtools stats' jobs in batch list, should be just 1"
        assert (
            batch.select_jobs('samtools stats')[0].name == j.name
        ), "Job name does not match job in batch"
        assert j is not None, "The samtools_stats function did not create a job."

        assert j._image == image_path('samtools')  # Ensure the correct image is used

    def test_sample_command(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        samtools_jobs = batch.select_jobs('samtools stats')
        # ---- Assertions
        cmd = get_command_str(samtools_jobs[0])
        print()

    def test_check_batch_command_executed(self, tmp_path: Path):
        pass

    def test_check_output_files_exist(self, tmp_path: Path):
        pass
