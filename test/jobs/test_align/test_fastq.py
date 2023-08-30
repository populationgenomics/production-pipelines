import re
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.jobs.align import Aligner, MarkDupTool, align

from ... import set_config
from ...factories.alignment_input import create_fastq_pairs_input
from ...factories.batch import create_local_batch
from ...factories.config import PipelineConfig
from ...factories.sequencing_group import create_sequencing_group

from ..helpers import get_command_str

from .shared import select_jobs, default_config


ALIGNERS = [Aligner.BWA, Aligner.BWAMEM2, Aligner.DRAGMAP]


class TestFastqPairsInput:
    """
    Test the `align` function using a SequencingGroup with a FASTQ pair(s)
    alignment input.
    """

    @staticmethod
    def _setup(config: PipelineConfig, tmp_path: Path, num_pairs: int = 1):
        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)

        dataset_id = config.workflow.dataset
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(
                prefix='SAMPLE1', location=tmp_path, n=num_pairs
            ),
        )

        return config, batch, sg

    @pytest.mark.parametrize('n', [1, 10])
    @pytest.mark.parametrize('aligner', ALIGNERS)
    def test_creates_one_align_job_for_each_pair(
        self, n: int, aligner: Aligner, tmp_path: Path
    ):
        config = default_config()
        config, batch, sg = self._setup(config, tmp_path, num_pairs=n)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == n

    @pytest.mark.parametrize('n', [1, 2])
    @pytest.mark.parametrize('aligner', ALIGNERS)
    def test_creates_merge_job_if_pairs_greater_than_one(
        self, n: int, aligner: Aligner, tmp_path: Path
    ):
        config = default_config()
        config, batch, sg = self._setup(config, tmp_path, num_pairs=n)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        merge_jobs = select_jobs(jobs, 'merge')
        assert len(merge_jobs) == int(n > 1)

        for job in merge_jobs:
            cmd = get_command_str(job)
            pattern = r'samtools merge -@\d+ -'
            for i in range(n):
                pattern += fr' \${{BATCH_TMPDIR}}/.*_SAMPLE1_L{i+1}_R_1_2.*/sorted_bam'
            pattern += r' > \${BATCH_TMPDIR}/.*/sorted_bam'

            match = re.search(pattern, cmd)
            if n > 1:
                assert match is not None
            else:
                assert match is None

    @pytest.mark.parametrize('n', [1, 2])
    @pytest.mark.parametrize('aligner', ALIGNERS)
    def test_sorts_output_with_samtools_if_more_than_one_pair(
        self, tmp_path: Path, n: int, aligner: Aligner
    ):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=n)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == n

        for job in align_jobs:
            cmd = get_command_str(job)
            pattern = r'\| samtools sort .* -Obam -o \${BATCH_TMPDIR}/.*/sorted_bam'
            match = re.search(pattern, cmd)
            if n > 1:
                assert match is not None
            else:
                assert match is None

    def test_using_bwa_generates_alignment_correct_command(self, tmp_path: Path):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=1)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.BWA,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        ref_file = config.references['broad']['ref_fasta']  # type: ignore
        cmd = get_command_str(align_jobs[0])
        assert re.search('bwa mem -K 100000000', cmd)  # Base pair processing
        assert re.search('-Y', cmd)  # Use soft-clipping
        assert re.search(rf"-R '\@RG\\tID:{sg.id}\\tSM:{sg.id}'", cmd)
        assert re.search(fr'\${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)
        assert re.search(r'\$BATCH_TMPDIR/R1.fq.gz \$BATCH_TMPDIR/R2.fq.gz', cmd)
        assert re.search(
            r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam', cmd
        )

    def test_using_bwamem2_generates_alignment_correct_command(self, tmp_path: Path):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=1)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.BWAMEM2,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        ref_file = config.references['broad']['ref_fasta']  # type: ignore
        cmd = get_command_str(align_jobs[0])
        assert re.search('bwa-mem2 mem -K 100000000', cmd)  # Base pair processing
        assert re.search('-Y', cmd)  # Use soft-clipping
        assert re.search(rf"-R '\@RG\\tID:{sg.id}\\tSM:{sg.id}'", cmd)
        assert re.search(fr'\${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)
        assert re.search(r'\$BATCH_TMPDIR/R1.fq.gz \$BATCH_TMPDIR/R2.fq.gz', cmd)
        assert re.search(
            r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam', cmd
        )

    @pytest.mark.parametrize('aligner', [Aligner.BWA, Aligner.BWAMEM2])
    def test_bwa_and_bwa2_uses_workflow_ref_fasta_if_set_in_config(
        self, tmp_path: Path, aligner: Aligner
    ):
        config = default_config()
        config.workflow.ref_fasta = 'workflow_ref.fasta'
        _, batch, sg = self._setup(config, tmp_path, num_pairs=1)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        cmd = get_command_str(align_jobs[0])
        assert config.references['broad']['ref_fasta'] not in cmd  # type: ignore
        assert config.workflow.ref_fasta in cmd

    def test_using_dragmap_generates_alignment_correct_command(self, tmp_path: Path):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=1)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        cmd = get_command_str(align_jobs[0])
        assert re.search(r'dragen-os -r \${BATCH_TMPDIR}/inputs/\w+', cmd)
        assert re.search(r'-1 \$BATCH_TMPDIR/R1.fq.gz -2 \$BATCH_TMPDIR/R2.fq.gz', cmd)
        assert re.search(
            r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam', cmd
        )

    def test_dragmap_aligner_reads_dragmap_reference_resource(
        self, mocker: MockFixture, tmp_path: Path
    ):
        config = default_config()
        config, batch, sg = self._setup(config, tmp_path)

        spy = mocker.spy(batch, 'read_input_group')
        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        file_location = config.references['broad']['dragmap_prefix']  # type: ignore
        spy.assert_any_call(
            **{
                # Trailing slash already exists in file_location
                'hash_table_cfg_bin': f'{file_location}hash_table.cfg.bin',
                'hash_table_cmp': f'{file_location}hash_table.cmp',
                'reference_bin': f'{file_location}reference.bin',
            }
        )

    @pytest.mark.parametrize('aligner', ALIGNERS)
    def test_biobambam_uses_output_from_merge_job_when_multiple_align_jobs_are_created(
        self, tmp_path: Path, aligner: Aligner
    ):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=5)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.BIOBAMBAM,
        )

        markdup_jobs = select_jobs(jobs, 'merge')
        assert len(markdup_jobs) == 1

        cmd = get_command_str(markdup_jobs[0])
        ref = config.references['broad']['ref_fasta']  # type: ignore
        assert re.search('| bamsormadup inputformat=bam', cmd)
        assert re.search('SO=coordinate', cmd)
        assert re.search(r'M=\${BATCH_TMPDIR}/.*/markdup_metrics', cmd)
        assert re.search('outputformat=sam', cmd)

        # Test SAM output of markdup is converted to CRAM
        assert re.search(
            (
                fr'| samtools view @\d+ -T \${{BATCH_TMPDIR}}/inputs/\w/{ref} \\'
                + r'\n-Ocram -o \${BATCH_TMPDIR}/.*/output_cram.cram'
            ),
            cmd,
        )

        # Test CRAM conversion is indexed
        assert re.search(
            (
                r'samtools index -@\d+ \${BATCH_TMPDIR}/.*/output_cram.cram \\'
                + r'\n\${BATCH_TMPDIR}/.*/output_cram.cram.crai'
            ),
            cmd,
        )

    @pytest.mark.parametrize('aligner', ALIGNERS)
    def test_biobambam_uses_output_from_align_job_when_one_align_job_is_created(
        self, tmp_path: Path, aligner: Aligner
    ):
        config = default_config()
        _, batch, sg = self._setup(config, tmp_path, num_pairs=1)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.BIOBAMBAM,
        )

        markdup_jobs = select_jobs(jobs, 'align')
        assert len(markdup_jobs) == 1

        cmd = get_command_str(markdup_jobs[0])
        ref = config.references['broad']['ref_fasta']  # type: ignore
        assert re.search('| bamsormadup inputformat=sam', cmd)
        assert re.search('SO=coordinate', cmd)
        assert re.search(r'M=\${BATCH_TMPDIR}/.*/markdup_metrics', cmd)
        assert re.search('outputformat=sam', cmd)

        # Test SAM output of markdup is converted to CRAM
        assert re.search(
            (
                fr'| samtools view @\d+ -T \${{BATCH_TMPDIR}}/inputs/\w/{ref} \\'
                + r'\n-Ocram -o \${BATCH_TMPDIR}/.*/output_cram.cram'
            ),
            cmd,
        )

        # Test CRAM conversion is indexed
        assert re.search(
            (
                r'samtools index -@\d+ \${BATCH_TMPDIR}/.*/output_cram.cram \\'
                + r'\n\${BATCH_TMPDIR}/.*/output_cram.cram.crai'
            ),
            cmd,
        )
