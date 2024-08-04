"""
Test the `align` function using a SequencingGroup with a BAM/CRAM alignment input.
"""

import enum
import re
from pathlib import Path
from typing import Optional

import pytest
from pytest_mock import MockFixture

from cpg_workflows.filetypes import BamPath, CramPath
from cpg_workflows.jobs.align import (
    Aligner,
    MarkDupTool,
    MissingAlignmentInputException,
    align,
)

from ... import set_config
from ...factories.alignment_input import create_bam_input, create_cram_input
from ...factories.batch import create_local_batch
from ...factories.config import PipelineConfig, StorageConfig
from ...factories.sequencing_group import create_sequencing_group
from ..helpers import get_command_str
from .shared import default_config, select_jobs


def setup_test(
    config: PipelineConfig,
    tmp_path: Path,
    alignment_input: Optional[CramPath | BamPath | str] = None,
    reference_assembly: Optional[str] = 'input_ref.fa',
    index: bool = True,
    create_realignment_cram: bool = False,
):
    set_config(config, tmp_path / 'config.toml')
    batch = create_local_batch(tmp_path)

    if isinstance(alignment_input, str):
        if alignment_input == 'cram':
            alignment_input = create_cram_input(
                prefix='SAMPLE1',
                location=tmp_path,
                index=index,
                reference_assembly=reference_assembly,
            )
        elif alignment_input == 'bam':
            alignment_input = create_bam_input(
                prefix='SAMPLE1',
                location=tmp_path,
                index=index,
            )
        else:
            raise ValueError(f'alignment_input must be "cram" or "bam", not {alignment_input}')

    seq_group = create_sequencing_group(
        dataset=config.workflow.dataset,
        sequencing_type=config.workflow.sequencing_type,
        alignment_input=alignment_input,
    )

    if config.workflow.realign_from_cram_version and create_realignment_cram:
        # New cram will be built from dataset path without an index which needs to exist
        # for the realignment reference to be used. No sharding will be performed since
        # the index is not available.
        path = seq_group.dataset.prefix() / 'cram' / 'v1' / f'{seq_group.id}.cram'
        path.parent.mkdir(parents=True)
        path.touch()

    return batch, seq_group


class TestPreProcessing:
    @pytest.mark.parametrize('input_type', ['bam', 'cram'])
    def test_does_not_shard_if_sequencing_exome(self, tmp_path: Path, input_type: str):
        """
        Test that when exome sequencing is specified, the `align` function does not
        split the alignment into more than one job.
        """
        config = default_config()
        config.workflow.sequencing_type = 'exome'

        batch, sg = setup_test(config, tmp_path, alignment_input=input_type)

        jobs = align(b=batch, sequencing_group=sg)
        assert len(select_jobs(jobs, 'align')) == 1

    @pytest.mark.parametrize('input_type', ['bam', 'cram'])
    @pytest.mark.parametrize('index,expected_count', [(True, 10), (False, 1)])
    def test_does_not_shard_if_sequencing_genome_and_index_does_not_exist(
        self,
        tmp_path: Path,
        input_type: str,
        index: bool,
        expected_count: int,
    ):
        """
        Test that when genome sequencing is specified, the `align` function does not
        split the alignment into more than one job unless an index exists.
        """
        config = default_config()
        config.workflow.sequencing_type = 'genome'

        batch, sg = setup_test(config, tmp_path, alignment_input=input_type, index=index)

        jobs = align(b=batch, sequencing_group=sg)
        assert len(select_jobs(jobs, 'align')) == expected_count

    @pytest.mark.parametrize('input_type', ['bam', 'cram'])
    def test_extract_fastq_sets_correct_flags_for_input_type(self, tmp_path: Path, input_type: str):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input=input_type)

        align(b=batch, sequencing_group=sg)
        all_jobs = batch._jobs
        extract_fastq_jobs = [job for job in all_jobs if job.name == 'Extract fastq']
        assert len(extract_fastq_jobs) == 1
        extract_fastq_j = extract_fastq_jobs[0]
        assert extract_fastq_j is not None

        cmd = get_command_str(extract_fastq_j)
        file = 'SAMPLE1.bam' if input_type == 'bam' else 'SAMPLE1.cram'
        escaped_file = re.escape(file)
        samtools_command = fr'\${{BATCH_TMPDIR}}/inputs/\w+/{escaped_file}'
        assert re.search(samtools_command, cmd)

    def test_extract_fastq_uses_correct_reference_assembly_for_cram(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='cram', reference_assembly='non_default.fa')

        align(b=batch, sequencing_group=sg)
        all_jobs = batch._jobs
        extract_fastq_jobs = [job for job in all_jobs if job.name == 'Extract fastq']
        assert len(extract_fastq_jobs) == 1
        extract_fastq_j = extract_fastq_jobs[0]
        assert extract_fastq_j is not None

        cmd = get_command_str(extract_fastq_j)
        ref = re.escape('non_default.fa')
        assert re.search(fr'samtools collate --reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref}', cmd)

    def test_extract_fastq_does_not_use_reference_flag_for_bam(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='bam')

        align(b=batch, sequencing_group=sg)
        all_jobs = batch._jobs
        extract_fastq_jobs = [job for job in all_jobs if job.name == 'Extract fastq']
        assert len(extract_fastq_jobs) == 1
        extract_fastq_j = extract_fastq_jobs[0]
        assert extract_fastq_j is not None

        cmd = get_command_str(extract_fastq_j)
        assert 'samtools collate' in cmd
        assert '--reference' not in cmd

    @pytest.mark.parametrize(
        'realignment_config,create_realignment_cram,expected_cram,expected_ref',
        [
            (
                {
                    'realign_from_cram_version': 'v1',
                    'cram_version_reference': {'v1': 'realign.fa'},
                },
                True,
                'CPGAAAAAA.cram',
                'realign.fa',
            ),
            (
                {
                    'realign_from_cram_version': 'v1',
                    'cram_version_reference': {'v1': 'realign.fa'},
                },
                False,
                'SAMPLE1.cram',
                'cram_ref.fa',
            ),
            (
                {
                    'realign_from_cram_version': None,
                    'cram_version_reference': None,
                },
                False,
                'SAMPLE1.cram',
                'cram_ref.fa',
            ),
        ],
    )
    def test_extract_fastq_and_uses_correct_reference_and_cram_file_when_realigning(
        self,
        tmp_path: Path,
        realignment_config: dict[str, str],
        create_realignment_cram: bool,
        expected_cram: str,
        expected_ref: str,
    ):
        """
        Test that the `align` function uses bazam to extract reads from a CRAM input
        relative to the reference assembly specified in the CRAM input, or relative to
        the realignment reference assembly if specified in config and it exists on disk
        """
        config = default_config()
        config.set_storage(config.workflow.dataset, StorageConfig(default=tmp_path))
        config.workflow.realign_from_cram_version = realignment_config['realign_from_cram_version']
        config.workflow.cram_version_reference = realignment_config['cram_version_reference']  # type: ignore

        batch, sg = setup_test(
            config,
            tmp_path,
            alignment_input='cram',
            reference_assembly='cram_ref.fa',
            create_realignment_cram=create_realignment_cram,
            index=False,
        )

        jobs = align(b=batch, sequencing_group=sg)
        align_jobs = select_jobs(jobs, 'align')
        all_jobs = batch._jobs
        extract_fastq_jobs = [job for job in all_jobs if job.name == 'Extract fastq']
        assert len(extract_fastq_jobs) == 1
        assert len(align_jobs) == 1
        extract_fastq_j = extract_fastq_jobs[0]

        ref = re.escape(expected_ref)
        file = re.escape(expected_cram)
        cmd = get_command_str(extract_fastq_j)
        pattern = fr'samtools collate --reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref} -@15 -u -O \${{BATCH_TMPDIR}}/inputs/\w+/{file} \$BATCH_TMPDIR/collate'

        assert re.search(pattern, cmd)

    def test_error_if_no_reference_assembly_on_cram_input(self, tmp_path: Path):
        """
        Test that the `align` function throws an error if there is no reference assembly
        specified on the CRAM input.
        """
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='cram', reference_assembly=None)

        with pytest.raises(
            AssertionError,
            match=r'The reference input for the alignment input .* was not set',
        ):
            align(b=batch, sequencing_group=sg)

    @pytest.mark.parametrize(
        'alignment_input,expected_error',
        [
            (None, r'No alignment inputs found for sequencing group .*'),
            ('cram', r'Alignment inputs for sequencing group .* do not exist'),
        ],
    )
    def test_error_if_alignment_input_does_not_exist(self, tmp_path: Path, alignment_input: str, expected_error: str):
        """
        Test that the `align` function throws an error if alignment input is `None`
        or doesn't exist on disk.
        """
        config = default_config()
        config.workflow.check_inputs = True
        batch, sg = setup_test(config, tmp_path, alignment_input=alignment_input)

        with pytest.raises(MissingAlignmentInputException, match=expected_error):
            align(b=batch, sequencing_group=sg)


class TestDragmap:
    @pytest.mark.parametrize('input_type', ['bam', 'cram'])
    def test_using_dragmap_generates_alignment_correct_command(self, tmp_path: Path, input_type: str):
        """
        Tests that the dragmap aligner command is correctly configured.
        """
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input=input_type)

        jobs = align(
            b=batch,
            sequencing_group=sg,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 10

        pattern = (
            r'dragen-os -r \${BATCH_TMPDIR}/inputs/\w+ --interleaved=1 -b r1 \\'
            fr'\n--RGID {sg.id} --RGSM {sg.id}'
            r'.* \| samtools sort .* -Obam -o \${BATCH_TMPDIR}/.*/sorted_bam'
        )

        cmd = get_command_str(align_jobs[0])
        assert re.search(pattern, cmd, flags=re.DOTALL)

    def test_dragmap_aligner_reads_dragmap_reference_resource(self, mocker: MockFixture, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='cram')

        spy = mocker.spy(batch, 'read_input_group')
        jobs = align(b=batch, sequencing_group=sg, aligner=Aligner.DRAGMAP)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 10

        file_location = config.references['broad']['dragmap_prefix']
        spy.assert_any_call(
            **{
                # Trailing slash already exists in file_location
                'hash_table_cfg_bin': f'{file_location}hash_table.cfg.bin',
                'hash_table_cmp': f'{file_location}hash_table.cmp',
                'reference_bin': f'{file_location}reference.bin',
            },
        )
