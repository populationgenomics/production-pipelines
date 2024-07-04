"""
Test the `align` function using a SequencingGroup with a BAM/CRAM alignment input.
"""

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

    def test_bazam_flags_correctly_set_when_extracting_reads(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='bam')

        jobs = align(b=batch, sequencing_group=sg)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 10

        file = re.escape('SAMPLE1.bam')
        for i, job in enumerate(align_jobs):
            cmd = get_command_str(job)
            pattern = fr'bazam .* -bam \${{BATCH_TMPDIR}}/inputs/\w+/{file} -s {i+1},10 > r1'
            assert re.search(pattern, cmd)

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
    def test_bazam_and_uses_correct_reference_and_cram_file_when_realigning(
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
        assert len(align_jobs) == 1

        ref = re.escape(expected_ref)
        file = re.escape(expected_cram)
        cmd = get_command_str(align_jobs[0])
        pattern = (
            fr'bazam .* -Dsamjdk\.reference_fasta=\${{BATCH_TMPDIR}}/inputs/\w+/{ref}'
            fr'.* -bam \${{BATCH_TMPDIR}}/inputs/\w+/{file} > r1'
        )
        assert re.search(pattern, cmd)

    def test_indexes_bam_if_index_does_not_exist(self, tmp_path: Path):
        """
        Test that the `align` function sorts and indexes the BAM input when index is
        missing.
        """
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='bam', index=False)

        jobs = align(b=batch, sequencing_group=sg)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1

        cmd = get_command_str(align_jobs[0])
        file = re.escape('SAMPLE1.bam')

        # Test sorts and indexes alignment input,
        assert re.search(
            (fr'samtools sort \${{BATCH_TMPDIR}}/inputs/\w+/{file} .* > \$BATCH_TMPDIR/sorted\.bam'),
            cmd,
        )
        # Test override original input
        assert re.search(fr'mv \$BATCH_TMPDIR/sorted\.bam \${{BATCH_TMPDIR}}/inputs/\w+/{file}', cmd)
        # Test indexes sorted input
        assert re.search(
            (
                fr'alignment_path="\${{BATCH_TMPDIR}}/inputs/\w+/{file}"'
                + r'\nsamtools index -@\d+ \$alignment_path \${alignment_path%m}i'
            ),
            cmd,
        )

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
    def test_indexes_cram_if_index_does_not_exist_and_uses_correct_reference(
        self,
        tmp_path: Path,
        realignment_config: dict[str, str],
        create_realignment_cram: bool,
        expected_cram: str,
        expected_ref: str,
    ):
        """
        Test that the `align` function sorts and indexes the CRAM input when index is
        missing, relative to the reference assembly specified in the CRAM input, or
        relative to the realignment reference assembly if specified in config and it
        exists on disk
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
        assert len(align_jobs) == 1

        cmd = get_command_str(align_jobs[0])
        file = re.escape(expected_cram)
        ref = re.escape(expected_ref)
        # Test sorts and indexes alignment input relative to reference
        assert re.search(
            (
                fr'samtools sort --reference \${{BATCH_TMPDIR}}/inputs/\w+/{ref} '
                fr'\${{BATCH_TMPDIR}}/inputs/\w+/{file}'
                r' .* > \$BATCH_TMPDIR/sorted\.cram'
            ),
            cmd,
        )
        # Test override original input
        assert re.search(
            fr'mv \$BATCH_TMPDIR/sorted.cram \${{BATCH_TMPDIR}}/inputs/\w+/{file}',
            cmd,
        )
        # Test indexes sorted input
        assert re.search(
            (
                fr'alignment_path="\${{BATCH_TMPDIR}}/inputs/\w+/{file}"'
                + r'\nsamtools index -@\d+ \$alignment_path \${alignment_path%m}i'
            ),
            cmd,
        )

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


class TestBwaAndBw2:
    @pytest.mark.parametrize('input_type', ['bam', 'cram'])
    @pytest.mark.parametrize('alinger, expected_tool', [(Aligner.BWA, 'bwa'), (Aligner.BWAMEM2, 'bwa-mem2')])
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
        input_type: str,
        reference: dict[str, str | None],
        expected_reference: str,
        alinger: Aligner,
        expected_tool: str,
    ):
        """
        Test that the BWA aligner command is correctly configured. Test that for each
        combination of alignment input type and aligner, the correct reference from the
        config is used.
        """
        config = default_config()
        config.workflow.ref_fasta = reference['workflow']
        config.references['broad']['ref_fasta'] = reference['broad']
        batch, sg = setup_test(config, tmp_path, alignment_input=input_type)

        jobs = align(b=batch, sequencing_group=sg, aligner=alinger)
        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 10

        ref_file = re.escape(expected_reference)
        cmd = get_command_str(align_jobs[0])
        # Test using base pair chunks, smart pairing and soft-clipping
        assert re.search(fr'{expected_tool} mem -K 100000000.*-p .* -Y', cmd, re.DOTALL)
        assert re.search(rf"-R '\@RG\\tID:{sg.id}\\tSM:{sg.id}'", cmd)
        assert re.search(fr'\${{BATCH_TMPDIR}}/inputs/\w+/{ref_file} r1', cmd)
        assert re.search(r'\| samtools sort .* -Obam -o \${BATCH_TMPDIR}/.*/sorted_bam', cmd)


class TestPostProcess:
    def test_creates_merge_job_if_input_sharded(self, tmp_path: Path):
        config = default_config()
        batch, sg = setup_test(config, tmp_path, alignment_input='cram')

        jobs = align(b=batch, sequencing_group=sg)

        merge_jobs = select_jobs(jobs, 'merge')
        assert len(merge_jobs) == 1

        # If there is more than one job, the merge job should them merge all provided
        # as space-separated inputs.
        cmd = get_command_str(merge_jobs[0])
        n = len(select_jobs(jobs, 'align'))
        pattern = (
            r'samtools merge -@\d+ - '
            + ' '.join(n * [r'\${BATCH_TMPDIR}/[A-Za-z0-9_-]+/sorted_bam'])
            + r' > \${BATCH_TMPDIR}/.*/sorted_bam'
        )
        assert re.search(pattern, cmd)

    def test_final_output_is_sorted_if_input_not_sharded(self, tmp_path: Path):
        config = default_config()
        config.workflow.sequencing_type = 'exome'
        batch, sg = setup_test(config, tmp_path, alignment_input='cram')

        jobs = align(b=batch, sequencing_group=sg)

        align_jobs = select_jobs(jobs, 'align')
        assert len(align_jobs) == 1
        assert len(select_jobs(jobs, 'merge')) == 0

        cmd = get_command_str(align_jobs[0])
        assert re.search(r'\| samtools sort .* -Obam > \${BATCH_TMPDIR}/.*/sorted_bam', cmd)
