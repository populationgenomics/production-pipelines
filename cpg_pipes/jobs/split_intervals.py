"""
Create Hail Batch jobs to split genomics intervals for parallel variant calling.
"""

import logging

import hailtop.batch as hb
from cpg_utils.hail_batch import reference_path, image_path
from hailtop.batch.job import Job

from cpg_pipes import utils, Path
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.types import SequencingType
from cpg_pipes.hb.command import wrap_command

logger = logging.getLogger(__file__)


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    intervals_path: Path | str | None = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
    job_attrs: dict | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job, list[hb.Resource]]:
    """
    Add a job that split genome into intervals to parallelise variant calling.

    As input interval file, takes intervals_path if provided, otherwise checks refs
    for the intervals of provided sequencing_type.

    This job calls picard's IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredicted number of intervals. WDL can handle
    that, but Hail Batch is not dynamic and have to expect certain number of output
    files.
    """
    job_attrs = (job_attrs or {}) | dict(tool='picard_IntervalListTools')
    job_name = f'Make {scatter_count} intervals for {sequencing_type.value}'
    
    if output_prefix and utils.exists(output_prefix / '1.interval_list'):
        job_attrs['reuse'] = True
        return b.new_job(job_name, job_attrs), [
            b.read_input(str(output_prefix / f'{idx + 1}.interval_list'))
            for idx in range(scatter_count)
        ]

    j = b.new_job(job_name, job_attrs)

    if not intervals_path:
        # Did we cache split intervals for this sequencing_type?
        cache_bucket = (
            reference_path('intervals_prefix') / 
            sequencing_type.value / 
            f'{scatter_count}intervals'
        )
        if utils.exists(cache_bucket / '1.interval_list'):
            j.name += ' [use cached]'
            return j, [
                b.read_input(str(cache_bucket / f'{idx + 1}.interval_list'))
                for idx in range(scatter_count)
            ]
        # Taking intervals file for the sequencing_type.
        intervals_path = reference_path(
            f'broad/{sequencing_type.value}_calling_interval_lists',
        )

    j.image(image_path('picard'))
    STANDARD.set_resources(j, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        SequencingType.GENOME: 100000,
        SequencingType.EXOME: 0,
    }.get(sequencing_type, 0)

    cmd = f"""
    mkdir /io/batch/out
    
    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(intervals_path))} \
    OUTPUT=/io/batch/out
    ls /io/batch/out
    ls /io/batch/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln /io/batch/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(wrap_command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )
    return j, [j[f'{idx + 1}.interval_list'] for idx in range(scatter_count)]
