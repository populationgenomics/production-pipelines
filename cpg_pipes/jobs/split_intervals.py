"""
Create Hail Batch jobs to split genomics intervals for parallel variant calling.
"""

import logging

import hailtop.batch as hb

from cpg_pipes import Path
from cpg_pipes import images
from cpg_pipes.refdata import RefData
from cpg_pipes.hb.command import wrap_command

logger = logging.getLogger(__file__)


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    refs: RefData,
    out_bucket: Path | None = None,
) -> hb.ResourceGroup:
    """
    Add a job that split genome into intervals to parallelise GnarlyGenotyper

    Returns a ResourceGroup instead of a Job because if the intervals are already
    pre-computed, no need to submit a job.
    """
    if scatter_count in refs.precomputed_intervals:
        intervals_bucket = refs.precomputed_intervals[scatter_count]
        return b.read_input_group(**{
            f'interval_{idx}': f'{intervals_bucket}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        })

    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(images.GATK_IMAGE)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )
    unpadded_intervals = b.read_input(refs.unpadded_intervals)
    ref = refs.fasta_res_group(b)

    cmd = f"""
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
    -L {unpadded_intervals} \\
    -O {j.intervals} \\
    -scatter {scatter_count} \\
    -R {ref.base} \\
    -mode INTERVAL_SUBDIVISION
    """
    j.command(wrap_command(cmd))
    if out_bucket:
        b.write_output(j.intervals, str(out_bucket / f'{scatter_count}intervals'))
    return j.intervals
