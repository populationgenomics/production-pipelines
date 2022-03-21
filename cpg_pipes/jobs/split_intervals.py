"""
Create Hail Batch jobs to split genomics intervals for parallel variant calling.
"""

import logging

import hailtop.batch as hb

from cpg_pipes import Path
from cpg_pipes import images
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.types import SequencingType
from cpg_pipes.refdata import RefData
from cpg_pipes.hb.command import wrap_command

logger = logging.getLogger(__file__)


def get_intervals(
    b: hb.Batch,
    refs: RefData,
    sequencing_type: SequencingType,
    scatter_count: int,
    out_bucket: Path | None = None,
) -> hb.ResourceGroup:
    """
    Add a job that split genome into intervals to parallelise GnarlyGenotyper

    Returns a ResourceGroup instead of a Job because if the intervals are already
    pre-computed, no need to submit a job.
    
    This job calls picard's IntervalListTools to scatter the input interval list 
    into scatter_count sub-interval lists.
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image(images.GATK_IMAGE)
    STANDARD.request_resources(storage_gb=16, mem_gb=2)
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    intervals = refs.calling_interval_lists[sequencing_type]
    break_bands_at_multiples_of = {
        SequencingType.WGS: 100000,
        SequencingType.EXOME: 0,
    }[sequencing_type]

    cmd = f"""
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.

    java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(intervals)} \
    OUTPUT={j.intervals}    
    """
    j.command(wrap_command(cmd))
    if out_bucket:
        b.write_output(j.intervals, str(out_bucket / f'{scatter_count}intervals'))
    return j.intervals
