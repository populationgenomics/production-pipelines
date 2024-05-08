"""
Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import command
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.jobs.bam_to_cram import cram_to_bam
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import (
    SequencingGroup,
)


def count_res_group(b: hb.Batch) -> hb.ResourceGroup:
    """
    Define resource group for counting.
    """
    gtf_file = get_config()['references']['star'].get('gtf')
    gtf_file = to_path(gtf_file)
    g = {
        'gtf': str(gtf_file),
    }
    return b.read_input_group(**g)


class FeatureCounts:
    """
    Construct a featureCounts command for counting reads.
    """

    def __init__(
        self,
        input_bam: BamPath | str | Path,
        gtf_file: str | Path,
        output_path: str | Path,
        summary_path: str | Path,
        paired_end: bool = True,
        feature_type: str = 'exon',
        attribute: str = 'gene_id',
        strandness: str = 'reverse',  # one of: 'none', 'forward', 'reverse'
        multi_mapping: bool = False,
        min_quality: int | None = None,
        primary_only: bool = False,
        ignore_duplicates: bool = False,
        count_pairs: bool = True,
        both_ends_mapped: bool = True,
        both_ends_same_chr: bool = True,
        threads: int = 1,
    ) -> None:
        self.command = [
            'featureCounts',
            *('-t', feature_type),
            *('-g', attribute),
            *('-s', {'none': '0', 'forward': '1', 'reverse': '2'}[strandness]),
            *('-a', str(gtf_file)),
            *('-T', str(threads)),
        ]

        if paired_end:
            self.command.append('-p')

        if multi_mapping:
            self.command.append('-M')

        if min_quality is not None:
            self.command.extend(['-Q', str(min_quality)])

        if primary_only:
            self.command.append('--primary')

        if ignore_duplicates:
            self.command.append('--ignoreDup')

        if paired_end and count_pairs:
            self.command.append('--countReadPairs')

        if both_ends_mapped:
            self.command.append('-B')

        if not both_ends_same_chr:
            self.command.append('-C')

        self.tmp_output = '$BATCH_TMPDIR/count_out/count'
        self.tmp_output_summary = f'{self.tmp_output}.summary'

        self.command.extend(['-o', self.tmp_output, str(input_bam)])

        self.make_tmpdir_command = 'mkdir -p $BATCH_TMPDIR/count_out'

        self.finalise_outputs_command = (
            f'ln {self.tmp_output} {output_path} && ln {self.tmp_output_summary} {summary_path}'
        )

    def __str__(self):
        return ' && '.join(
            [
                self.make_tmpdir_command,
                ' '.join(self.command),
                self.finalise_outputs_command,
            ],
        )

    def __repr__(self):
        return self.__str__()


def count(
    b: hb.Batch,
    input_cram_or_bam: BamPath | CramPath,
    output_path: str | Path,
    summary_path: str | Path,
    cram_to_bam_path: Path | None = None,
    sample_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> list[Job]:
    """
    Count RNA seq reads mapping to genes and/or transcripts using featureCounts.
    """
    # Reuse existing output if possible
    if output_path and summary_path and can_reuse(output_path, overwrite) and can_reuse(summary_path, overwrite):
        return []

    jobs: list[Job] = []

    # Determine whether input is a BAM file or a CRAM file
    if isinstance(input_cram_or_bam, BamPath):
        # Localise input
        input_bam_reads = b.read_input_group(
            **{
                'bam': str(input_cram_or_bam.path),
                'bam.bai': str(input_cram_or_bam.index_path),
            },
        )
    elif isinstance(input_cram_or_bam, CramPath):
        # Localise input
        input_cram_reads = b.read_input_group(
            **{
                'cram': str(input_cram_or_bam.path),
                'cram.crai': str(input_cram_or_bam.index_path),
            },
        )
        # Convert CRAM to BAM
        j, input_bam_reads = cram_to_bam(
            b=b,
            input_cram=input_cram_reads,
            output_bam=cram_to_bam_path,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        if j and isinstance(j, Job):
            jobs.append(j)
    else:
        raise ValueError(f'Invalid alignment input: "{str(input_cram_or_bam)}", expected BAM or CRAM file.')

    assert isinstance(input_bam_reads, hb.ResourceGroup)

    counting_reference = count_res_group(b)

    # Create job
    job_name = f'count_{sample_name}' if sample_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='featureCounts')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('subread'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,  # TODO: make configurable
    )

    # Declare output resource group
    j.declare_resource_group(
        count_output={
            'count': '{root}.count',
            'count.summary': '{root}.count.summary',
        },
    )

    # Create counting command
    fc = FeatureCounts(
        input_bam=input_bam_reads.bam,
        gtf_file=counting_reference.gtf,
        output_path=j.count_output['count'],
        summary_path=j.count_output['count.summary'],
        paired_end=True,
        feature_type='exon',
        attribute='gene_id',
        strandness='reverse',  # TODO: confirm this
        multi_mapping=False,  # TODO: determine default value
        min_quality=None,  # TODO: determine default value
        primary_only=True,  # TODO: determine default value
        ignore_duplicates=False,  # TODO: determine default value
        count_pairs=True,  # TODO: determine default value
        both_ends_mapped=True,  # TODO: determine default value
        both_ends_same_chr=True,  # TODO: determine default value
        threads=res.get_nthreads(),
    )
    cmd = str(fc)

    # Add command to job
    j.command(command(cmd, monitor_space=True))

    jobs.append(j)

    # Write output to file
    if output_path:
        b.write_output(j.count_output['count'], str(output_path))
    if summary_path:
        b.write_output(j.count_output['count.summary'], str(summary_path))

    return jobs
