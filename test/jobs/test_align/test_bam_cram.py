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

from .shared import select_jobs, default_config


ALIGNERS = [Aligner.BWA, Aligner.BWAMEM2, Aligner.DRAGMAP]

pytestmark = pytest.mark.parametrize('aligner', ALIGNERS)


def _setup(
    config: PipelineConfig,
    tmp_path: Path,
    alignment_input: Optional[CramPath | BamPath] = None,
):
    set_config(config, tmp_path / 'config.toml')
    batch = create_local_batch(tmp_path)

    if alignment_input is None:
        alignment_input = create_bam_input(prefix='SAMPLE1', location=tmp_path)

    sg = create_sequencing_group(
        dataset=config.workflow.dataset,
        sequencing_type=config.workflow.sequencing_type,
        alignment_input=alignment_input,
    )

    return config, batch, sg


def test_does_not_shard_and_merge_realignment_if_sequencing_exome(
    tmp_path: Path, aligner: Aligner
):
    config = default_config()
    config.workflow.sequencing_type = 'exome'

    config, batch, sg = _setup(config, tmp_path)

    jobs = align(
        b=batch,
        sequencing_group=sg,
        aligner=aligner,
        markdup_tool=MarkDupTool.NO_MARKDUP,
    )

    assert len(select_jobs(jobs, 'align')) == 1
    assert len(select_jobs(jobs, 'merge')) == 0


def test_shards_and_merges_realignment_if_sequencing_genome(
    tmp_path: Path, aligner: Aligner
):
    """
    Test that the `align` function shards the alignment if the sequencing type is
    exome.
    """
    config = default_config()
    config.workflow.sequencing_type = 'genome'

    config, batch, sg = _setup(config, tmp_path)

    jobs = align(
        b=batch,
        sequencing_group=sg,
        aligner=aligner,
        markdup_tool=MarkDupTool.NO_MARKDUP,
    )

    align_jobs = select_jobs(jobs, 'align')
    assert len(align_jobs) == 10
    assert len(select_jobs(jobs, 'merge')) == 1

    for i, job in enumerate(align_jobs):
        cmd = get_command_str(job)
        assert re.search(
            fr'bazam .* -bam \${{BATCH_TMPDIR}}/inputs/\w+/SAMPLE1.bam -s {i+1},10',
            cmd,
        )


def test_does_not_shard_bam_and_indexes_bam_if_index_does_not_exist(
    tmp_path: Path, aligner: Aligner
):
    """
    Test that the `align` function does not shard a BAM input if it's index is missing,
    and instead sorts and indexes the BAM input.
    """
    config = default_config()
    config.workflow.sequencing_type = 'genome'

    config, batch, sg = _setup(
        config,
        tmp_path,
        alignment_input=create_bam_input(
            prefix='SAMPLE1', location=tmp_path, index=False
        ),
    )

    jobs = align(
        b=batch,
        sequencing_group=sg,
        aligner=aligner,
        markdup_tool=MarkDupTool.NO_MARKDUP,
    )

    align_jobs = select_jobs(jobs, 'align')
    assert len(align_jobs) == 1

    # Test sorts and indexes alignment input, overriding original.
    cmd = get_command_str(align_jobs[0])
    file = 'SAMPLE1.bam'
    assert re.search(fr'samtools sort \${{BATCH_TMPDIR}}/inputs/\w+/{file}', cmd)
    assert re.search(
        fr'mv \$BATCH_TMPDIR/sorted.bam \${{BATCH_TMPDIR}}/inputs/\w+/{file}', cmd
    )
    assert re.search(
        (
            fr'alignment_path="\${{BATCH_TMPDIR}}/inputs/\w+/{file}"\n'
            + r'samtools index -@\d+ \$alignment_path \${alignment_path%m}i'
        ),
        cmd,
    )


def test_does_not_shard_cram_and_indexes_cram_if_index_does_not_exist(
    tmp_path: Path, aligner: Aligner
):
    """
    Test that the `align` function does not shard a CRAM input if it's index is missing,
    and instead sorts and indexes the BAM input.
    """
    config = default_config()
    config.workflow.sequencing_type = 'genome'

    config, batch, sg = _setup(
        config,
        tmp_path,
        alignment_input=create_cram_input(
            prefix='SAMPLE1',
            location=tmp_path,
            reference_assembly='cram_ref.fa',
            index=False,
        ),
    )

    jobs = align(
        b=batch,
        sequencing_group=sg,
        aligner=aligner,
        markdup_tool=MarkDupTool.NO_MARKDUP,
    )

    align_jobs = select_jobs(jobs, 'align')
    assert len(align_jobs) == 1

    # Test sorts and indexes alignment input, overriding original.
    cmd = get_command_str(align_jobs[0])
    file = 'SAMPLE1.cram'
    assert re.search(
        (
            fr'samtools sort --reference \${{BATCH_TMPDIR}}/inputs/\w+/cram_ref.fa '
            + fr'\${{BATCH_TMPDIR}}/inputs/\w+/{file}'
        ),
        cmd,
    )
    assert re.search(
        fr'mv \$BATCH_TMPDIR/sorted.cram \${{BATCH_TMPDIR}}/inputs/\w+/{file}',
        cmd,
    )
    assert re.search(
        (
            fr'alignment_path="\${{BATCH_TMPDIR}}/inputs/\w+/{file}"\n'
            + r'samtools index -@\d+ \$alignment_path \${alignment_path%m}i'
        ),
        cmd,
    )


@pytest.mark.parametrize('exists', [True, False])
def test_use_reference_assembly_from_config_when_realigning_cram_version_if_exists(
    tmp_path: Path, aligner: Aligner, exists: bool
):
    """
    Test that the `align` function re-aligns a CRAM file using the re-alignment
    specified in the configration. If this file does not exist, then default to the
    reference assembly specified in the CRAM input.
    """
    config = default_config()
    config.set_storage(config.workflow.dataset, StorageConfig(default=tmp_path))
    config.workflow.realign_from_cram_version = 'v1'
    config.workflow.cram_version_reference = {'v1': 'realign.fa'}

    config, batch, sg = _setup(
        config,
        tmp_path,
        alignment_input=create_cram_input(
            prefix='SAMPLE1',
            location=tmp_path,
            reference_assembly='cram_ref.fa',
            index=False,
        ),
    )

    ref = sg.dataset.prefix() / 'cram' / 'v1' / f'{sg.id}.cram'
    ref.parent.mkdir(parents=True)
    if exists:
        ref.touch()

    jobs = align(
        b=batch,
        sequencing_group=sg,
        aligner=aligner,
        markdup_tool=MarkDupTool.NO_MARKDUP,
    )

    align_jobs = select_jobs(jobs, 'align')
    assert len(align_jobs) == 1

    cmd = get_command_str(align_jobs[0])
    if exists:
        assert re.search(
            fr'samtools sort --reference \${{BATCH_TMPDIR}}/inputs/\w+/realign.fa',
            cmd,
        )
    else:
        assert re.search(
            fr'samtools sort --reference \${{BATCH_TMPDIR}}/inputs/\w+/cram_ref.fa',
            cmd,
        )


def test_error_if_no_reference_assembly_on_cram_input(tmp_path: Path, aligner: Aligner):
    """
    Test that the `align` function throws an error if there is no reference assembly
    specified in the CRAM input.
    """
    config = default_config()

    config, batch, sg = _setup(
        config,
        tmp_path,
        alignment_input=create_cram_input(
            prefix='SAMPLE1', location=tmp_path, reference_assembly=None, index=False
        ),
    )

    with pytest.raises(
        AssertionError,
        match=r'The reference input for the alignment input .* was not set',
    ):
        align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )


def test_error_if_alignment_input_does_not_exist(tmp_path: Path, aligner: Aligner):
    """
    Test that the `align` function throws an error if there is no reference assembly
    specified in the CRAM input.
    """
    config = default_config()
    config.workflow.check_inputs = True

    config, batch, sg = _setup(config, tmp_path)

    with pytest.raises(
        MissingAlignmentInputException,
        match=r'Alignment inputs for sequencing group .* do not exist',
    ):
        align(
            b=batch,
            sequencing_group=sg,
            aligner=aligner,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )


# def test_generates_bwa
# def test_generates_bwamem2
# def test_uses_workflow_reference
# def test_error_no_broad_or_workflow_reference

# def test_generates_dragmap
# def test_error_no_dragmap_reference

# def test_generates_markdup
