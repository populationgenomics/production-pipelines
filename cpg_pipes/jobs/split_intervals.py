"""
Create Hail Batch jobs to split genomics intervals for parallel variant calling.
"""

import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

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
    job_attrs: dict | None = None,
) -> Job:
    """
    Add a job that split genome into intervals to parallelise variant calling.

    This job calls picard's IntervalListTools to scatter the input interval list 
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl
    
    Note that we use the mode INTERVAL_SUBDIVISION instead of 
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than 
    INTERVAL_SUBDIVISION produce an unpredicted number of intervals. WDL can handle
    that, but Hail Batch is not dynamic and have to expect certain number of output
    files.
    """
    j = b.new_job(f'Make {scatter_count} intervals', job_attrs)
    j.image(images.SAMTOOLS_PICARD_IMAGE)
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
    mkdir /io/batch/out
    
    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(intervals))} \
    OUTPUT=/io/batch/out
    ls /io/batch/out
    ls /io/batch/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
    ln /io/batch/out/{name}/scattered.interval_list {j[f'{idx}.interval_list']}
    """
    
    j.command(wrap_command(cmd))
    if out_bucket:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx}.interval_list'], 
                str(out_bucket / f'{scatter_count}intervals' / f'{idx}.interval_list')
            )
    return j
