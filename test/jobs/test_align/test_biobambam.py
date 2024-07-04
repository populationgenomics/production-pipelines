"""
TODO: This test should be moved to a markdup test module. Picard is defined in it's own
module, however biobambam2 is constructed in the align job module so it can use piped
output from an alignment or samtools merge command.
"""

import re
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs.align import MarkDupTool, align

from ...factories.config import PipelineConfig
from ..helpers import get_command_str
from .shared import default_config, select_jobs
from .test_bam_cram import setup_test as setup_bam_cram_test
from .test_fastq import setup_test as setup_fastq_test


def setup_test(config: PipelineConfig, tmp_path: Path, input_type: str, **kwargs):
    if input_type == 'fastq':
        # Remove kwargs not in setup_fastq_test
        kwargs.pop('alignment_input', None)
        return setup_fastq_test(config, tmp_path, **kwargs)
    else:
        # Remove kwargs not in setup_bam_cram_test
        kwargs.pop('num_pairs', None)
        return setup_bam_cram_test(config, tmp_path, alignment_input=input_type, **kwargs)


@pytest.mark.parametrize('input_type', ['fastq', 'bam', 'cram'])
def test_uses_output_from_merge_job_when_multiple_align_jobs_are_created(tmp_path: Path, input_type: str):
    config = default_config()
    batch, sg = setup_test(config, tmp_path, input_type=input_type, num_pairs=5)

    jobs = align(
        b=batch,
        sequencing_group=sg,
        markdup_tool=MarkDupTool.BIOBAMBAM,
    )

    markdup_jobs = select_jobs(jobs, 'merge')
    assert len(markdup_jobs) == 1

    cmd = get_command_str(markdup_jobs[0])
    ref = config.references['broad']['ref_fasta']
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


@pytest.mark.parametrize('input_type', ['fastq', 'bam', 'cram'])
def test_uses_output_from_align_job_when_one_align_job_is_created(tmp_path: Path, input_type: str):
    config = default_config()
    config.workflow.sequencing_type = 'exome'
    batch, sg = setup_test(config, tmp_path, input_type=input_type, num_pairs=1)

    jobs = align(
        b=batch,
        sequencing_group=sg,
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


@pytest.mark.parametrize('input_type', ['fastq', 'bam', 'cram'])
@pytest.mark.parametrize(
    'tool,expected_resource',
    [(MarkDupTool.NO_MARKDUP, 'sorted_bam'), (MarkDupTool.BIOBAMBAM, 'output_cram')],
)
def test_writes_correct_resource_depending_on_markdup_tool(
    mocker: MockFixture,
    tmp_path: Path,
    input_type: str,
    tool: MarkDupTool,
    expected_resource: str,
):
    config = default_config()
    batch, sg = setup_test(config, tmp_path, input_type=input_type, num_pairs=1)
    out = CramPath(path=tmp_path / 'out')

    spy = mocker.spy(batch, 'write_output')

    jobs = align(b=batch, sequencing_group=sg, output_path=out, markdup_tool=tool)

    final_job = jobs[-1]
    spy.assert_called_with(getattr(final_job, expected_resource), str(out.path.with_suffix('')))


@pytest.mark.parametrize('input_type', ['fastq', 'bam', 'cram'])
def test_writes_markdup_metrics_resrouce(mocker: MockFixture, tmp_path: Path, input_type: str):
    config = default_config()
    batch, sg = setup_test(config, tmp_path, input_type=input_type, num_pairs=1)
    cram_out = CramPath(path=tmp_path / 'out')
    metrics_out = Path(tmp_path / 'out')

    spy = mocker.spy(batch, 'write_output')
    jobs = align(
        b=batch,
        sequencing_group=sg,
        output_path=cram_out,
        out_markdup_metrics_path=metrics_out,
    )

    final_job = jobs[-1]
    spy.assert_called_with(final_job.markdup_metrics, str(metrics_out))
