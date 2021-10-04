import logging
from os.path import join
from typing import Optional

import hailtop.batch as hb

from cpg_production_pipelines import resources
from cpg_production_pipelines.jobs import wrap_command

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def intervals(
    b: hb.Batch,
    scatter_count: int,
    ref_fasta: str,
    out_bucket: Optional[str] = None,
) -> hb.ResourceGroup:
    """
    Add a job that split genome into intervals to parallelise GnarlyGenotyper

    Returns a ResourceGroup instead of a Job because if the intervals are already
    pre-computed, no need to submit a job.
    """
    if scatter_count in resources.PRECOMPUTED_INTERVALS:
        intervals_bucket = resources.PRECOMPUTED_INTERVALS[scatter_count]
        return b.read_input_group(
            intervals={
                f'interval_{idx}': f'{intervals_bucket}/{str(idx).zfill(4)}-scattered.interval_list'
                for idx in range(scatter_count)
            }
        )

    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(resources.GATK_IMAGE)
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    cmd = f"""
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
    -L {resources.UNPADDED_INTERVALS} \\
    -O {j.intervals} \\
    -scatter {scatter_count} \\
    -R {ref_fasta} \\
    -mode INTERVAL_SUBDIVISION
    """
    output_bucket_path_to_check = None
    if out_bucket:
        out_bucket = join(out_bucket, f'{scatter_count}intervals')
        # Checking the first and the last interval:
        output_bucket_path_to_check = [
            join(out_bucket, f'0000-scattered.interval_list'),
            join(
                out_bucket, f'{str(scatter_count-1).zfill(4)}-scattered.interval_list'
            ),
        ]
    j.command(wrap_command(cmd, output_bucket_path_to_check))
    if out_bucket:
        b.write_output(j.intervals, out_bucket)
    return j.intervals
