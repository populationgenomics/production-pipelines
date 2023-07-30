import re
from pathlib import Path

from cpg_workflows.jobs.align import Aligner, MarkDupTool, align

from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='align-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'dragmap': 'dragmap_image:1.3.0',
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


class TestSingleFastqPairAlignment:
    def test_creates_one_align_job(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
        )

        # ---- The job that we want to test
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        assert len(align_jobs) == 1

    def test_sorts_output_with_bamtools(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config.as_dict(), tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
        )

        # ---- The job that we want to test
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        cmd = get_command_str(align_jobs[0])
        assert re.search(r'\| samtools sort .* -Obam', cmd)
