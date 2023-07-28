import re
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.jobs.align import Aligner, MarkDupTool, align

from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from ..factories.types import SequencingType
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
            'dragmap': 'dragmap:latest',
            'bwa': 'bwa:latest',
            'bwamem2': 'bwamem2:latest',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
            }
        },
    )


# def test_error_no_dragmap_reference_when_using_dragmap_aligner


class TestSingleFastqPairAlignment:
    def test_creates_one_align_job(self, tmp_path: Path):
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
        )

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

    def test_dragen_aligner_sets_reads_pairs_in_correct_order(self, tmp_path: Path):
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        fastq = create_fastq_pairs_input(location=tmp_path, n=1)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=fastq,
        )

        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        expected_cmd = (
            r"dragen-os -r \${BATCH_TMPDIR}/inputs/\w+ "
            + r"-1 \$BATCH_TMPDIR/R1.fq.gz -2 \$BATCH_TMPDIR/R2.fq.gz"
        )
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        cmd = cmd = get_command_str(align_jobs[0])
        assert re.search(expected_cmd, cmd)

    def test_dragen_aligner_reads_dragmap_reference_resources(
        self, mocker: MockFixture, tmp_path: Path
    ):
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        fastq = create_fastq_pairs_input(location=tmp_path, n=1)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=fastq,
        )

        spy = mocker.spy(batch, 'read_input_group')
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        file_location = config.references['broad']['dragmap_prefix']  # type: ignore
        spy.assert_any_call(
            **{
                # Trailing slash already exists in file_location
                'hash_table_cfg_bin': f"{file_location}hash_table.cfg.bin",
                'hash_table_cmp': f"{file_location}hash_table.cmp",
                'reference_bin': f"{file_location}reference.bin",
            }
        )

    def test_using_bwa_generates_correct_command(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        fastq = create_fastq_pairs_input(location=tmp_path, n=1)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=fastq,
        )

        # ---- Run the job
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.BWA,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        ref_file = config.dig('references', 'broad', 'ref_fasta')
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        cmd = cmd = get_command_str(align_jobs[0])
        assert re.search('bwa mem -K 100000000', cmd)  # Base processing
        assert re.search('-Y', cmd)  # Use soft-clipping
        assert re.search(rf"-R '\@RG\\tID:{sg.id}\\tSM:{sg.id}'", cmd)
        assert re.search(fr"\${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}", cmd)
        assert re.search(r'\$BATCH_TMPDIR/R1.fq.gz \$BATCH_TMPDIR/R2.fq.gz', cmd)

    @pytest.mark.parametrize('seq_type', ['genome', 'exome'])
    @pytest.mark.parametrize('aligner', [Aligner.BWA, Aligner.BWAMEM2, Aligner.DRAGMAP])
    def test_sorts_output_with_bamtools(
        self, tmp_path: Path, aligner: Aligner, seq_type: SequencingType
    ):
        # ---- Test setup
        config = default_config()
        config.workflow.sequencing_type = seq_type
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
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        cmd = get_command_str(align_jobs[0])
        assert re.search(
            r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam', cmd
        )
