"""
Test the `align` function using a SequencingGroup with a FASTQ pair(s)
alignment input.
"""

import re
from calendar import c
from pathlib import Path
from typing import Optional

import pytest
from pytest_mock import MockFixture

from cpg_utils.config import ConfigError
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs.align import (
    Aligner,
    MarkDupTool,
    MissingAlignmentInputException,
    align,
)

from ... import set_config
from ...factories.alignment_input import create_fastq_pairs_input
from ...factories.batch import create_local_batch
from ...factories.config import PipelineConfig
from ...factories.sequencing_group import create_sequencing_group
from ..helpers import get_command_str
from .shared import default_config, select_jobs


def setup_test(
    config: PipelineConfig,
    tmp_path: Path,
    num_pairs: Optional[int] = 1,
    create=False,
):
    set_config(config, tmp_path / 'config.toml')
    batch = create_local_batch(tmp_path)

    dataset_id = config.workflow.dataset
    alignment_input = None
    if num_pairs:
        alignment_input = create_fastq_pairs_input(
            prefix='SAMPLE1',
            location=tmp_path,
            n=num_pairs,
            create=create,
        )

    sg = create_sequencing_group(
        dataset=dataset_id,
        sequencing_type=config.workflow.sequencing_type,
        alignment_input=alignment_input,
    )

    return batch, sg


class TestPreProcess:
    @pytest.mark.parametrize('n', [1, 2])
    def test_creates_one_align_job_for_each_fq_pair(self, tmp_path: Path, n: int):
        """
        Test that the `align` function creates one align job for each pair of FASTQ
        files.
        """
        config = default_config()
        batch, sg = setup_test(config, tmp_path, num_pairs=n)

        jobs = align(b=batch, sequencing_group=sg)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == n

    @pytest.mark.parametrize(
        'n_pairs,expected_error',
        [
            (None, r'No alignment inputs found for sequencing group .*'),
            (2, r'Alignment inputs for sequencing group .* do not exist'),
        ],
    )
    def test_error_if_alignment_input_is_missing(
        self,
        tmp_path: Path,
        n_pairs: Optional[int],
        expected_error: str,
    ):
        """
        Test that the `align` function throws an error if the alignment input is
        `None` or the file doesn't exist.
        """
        config = default_config()
        config.workflow.check_inputs = True

        batch, sg = setup_test(config, tmp_path, num_pairs=n_pairs)

        with pytest.raises(MissingAlignmentInputException, match=expected_error):
            align(b=batch, sequencing_group=sg)

    def test_error_if_no_workflow_of_broad_reference_set(self, tmp_path: Path):
        config = default_config()
        config.workflow.ref_fasta = None
        config.references['broad']['ref_fasta'] = None

        batch, sg = setup_test(config, tmp_path)

        with pytest.raises(ConfigError, match=r'Key "ref_fasta" not found'):
            align(
                b=batch,
                sequencing_group=sg,
                markdup_tool=MarkDupTool.NO_MARKDUP,
            )


class TestDragmap:
    def test_using_dragmap_generates_alignment_correct_command(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, num_pairs=2)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 2

        pattern = (
            r'dragen-os -r \${BATCH_TMPDIR}/inputs/\w+'
            r' -1 \$BATCH_TMPDIR/R1\.fq\.gz -2 \$BATCH_TMPDIR/R2\.fq\.gz'
            fr'.*--RGID {sg.id} --RGSM {sg.id}'
            r'.* \| samtools sort .* -Obam -o \${BATCH_TMPDIR}/.*/sorted_bam'
        )

        for job in align_jobs:
            cmd = get_command_str(job)
            assert re.search(pattern, cmd, flags=re.DOTALL)

    def test_dragmap_aligner_reads_dragmap_reference_resource(
        self,
        mocker: MockFixture,
        tmp_path: Path,
    ):
        config = default_config()
        batch, sg = setup_test(config, tmp_path)

        spy = mocker.spy(batch, 'read_input_group')
        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        file_location = config.references['broad']['dragmap_prefix']
        spy.assert_any_call(
            **{
                # Trailing slash already exists in file_location
                'hash_table_cfg_bin': f'{file_location}hash_table.cfg.bin',
                'hash_table_cmp': f'{file_location}hash_table.cmp',
                'reference_bin': f'{file_location}reference.bin',
            },
        )


class TestBwaAndBwa2:
    @pytest.mark.parametrize(
        'alinger,expected_tool',
        [(Aligner.BWA, 'bwa'), (Aligner.BWAMEM2, 'bwa-mem2')],
    )
    @pytest.mark.parametrize(
        'reference,expected_reference',
        [
            ({'workflow': 'workflow.fa', 'broad': 'broad.fa'}, 'workflow.fa'),
            ({'workflow': 'workflow.fa', 'broad': None}, 'workflow.fa'),
            ({'workflow': None, 'broad': 'broad.fa'}, 'broad.fa'),
        ],
    )
    def test_using_bwa_generates_alignment_correct_command(
        self,
        tmp_path: Path,
        alinger: Aligner,
        expected_tool: str,
        reference: dict[str, str | None],
        expected_reference: str,
    ):
        config = default_config()
        config.workflow.ref_fasta = reference['workflow']
        config.references['broad']['ref_fasta'] = reference['broad']
        batch, sg = setup_test(config, tmp_path, num_pairs=2)

        jobs = align(b=batch, sequencing_group=sg, aligner=alinger)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 2

        ref_file = re.escape(expected_reference)
        cmd = get_command_str(align_jobs[0])
        # Test using base pair chunks, smart pairing and soft-clipping
        assert re.search(fr'{expected_tool} mem -K 100000000.*-Y', cmd, re.DOTALL)
        assert re.search(rf"-R '\@RG\\tID:{sg.id}\\tSM:{sg.id}'", cmd)
        assert re.search(fr'\${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)
        assert re.search(r'\$BATCH_TMPDIR/R1.fq.gz \$BATCH_TMPDIR/R2.fq.gz', cmd)
        assert re.search(
            r'\| samtools sort .* -Obam -o \${BATCH_TMPDIR}/.*/sorted_bam',
            cmd,
        )


class TestPostProcess:
    def test_no_merge_job_create_for_single_fq_pair(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, num_pairs=1)

        jobs = align(b=batch, sequencing_group=sg)

        align_jobs = select_jobs(jobs, 'merge')
        assert len(align_jobs) == 0

    def test_creates_merge_job_if_fq_pairs_greater_than_one(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, num_pairs=2)

        jobs = align(b=batch, sequencing_group=sg)

        merge_jobs = select_jobs(jobs, 'merge')
        assert len(merge_jobs) == 1

        # If there is more than one pair, the merge job should them merge all provided
        # as space-separated inputs.
        cmd = get_command_str(merge_jobs[0])
        pattern = (
            r'samtools merge -@\d+ -'
            r' \${BATCH_TMPDIR}/.*_SAMPLE1_L1_R_1_2.*/sorted_bam'
            r' \${BATCH_TMPDIR}/.*_SAMPLE1_L2_R_1_2.*/sorted_bam'
            r' > \${BATCH_TMPDIR}/.*/sorted_bam'
        )
        assert re.search(pattern, cmd)

    def test_final_output_is_sorted_if_one_fastq_pair_is_provided(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, num_pairs=1)

        jobs = align(b=batch, sequencing_group=sg)

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        cmd = get_command_str(align_jobs[0])
        assert re.search(
            r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam',
            cmd,
        )
