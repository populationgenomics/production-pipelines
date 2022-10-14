"""
Utility jobs to parallelize VCF processing.
"""

from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, command
from cpg_utils.hail_batch import reference_path
from cpg_utils.workflows.resources import STANDARD
from cpg_utils.workflows.utils import can_reuse
from cpg_utils.workflows.utils import exists


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    source_intervals_path: Path | None = None,
    job_attrs: dict[str, str] | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job | None, list[hb.ResourceFile]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param source_intervals_path: path to source intervals to split. Would check for
        config if not provided.
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredictable number of intervals. WDL can
    handle that, but Hail Batch is not dynamic and expects a certain number
    of output files.
    """
    assert scatter_count > 0, scatter_count
    sequencing_type = get_config()['workflow']['sequencing_type']
    source_intervals_path = source_intervals_path or reference_path(
        f'broad/{sequencing_type}_calling_interval_lists'
    )

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(str(source_intervals_path))]

    if output_prefix and exists(output_prefix / '1.interval_list'):
        return None, [
            b.read_input(str(output_prefix / f'{idx + 1}.interval_list'))
            for idx in range(scatter_count)
        ]

    if not source_intervals_path and exists(
        (
            existing_split_intervals_prefix := (
                reference_path('intervals_prefix')
                / sequencing_type
                / f'{scatter_count}intervals'
            )
        )
        / '1.interval_list'
    ):
        # We already have split intervals for this sequencing_type:
        return None, [
            b.read_input(
                str(existing_split_intervals_prefix / f'{idx + 1}.interval_list')
            )
            for idx in range(scatter_count)
        ]

    j = b.new_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | dict(tool='picard_IntervalListTools'),
    )
    j.image(image_path('picard'))
    STANDARD.set_resources(j, storage_gb=16, mem_gb=2)

    break_bands_at_multiples_of = {
        'genome': 100000,
        'exome': 0,
    }.get(sequencing_type, 0)

    cmd = f"""
    mkdir $BATCH_TMPDIR/out
    
    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(source_intervals_path))} \
    OUTPUT=$BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln $BATCH_TMPDIR/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )

    intervals: list[hb.ResourceFile] = []
    for idx in range(scatter_count):
        interval = j[f'{idx + 1}.interval_list']
        assert isinstance(interval, hb.ResourceFile)
        intervals.append(interval)
    return j, intervals


def subset_vcf(
    b: hb.Batch,
    vcf: hb.ResourceGroup,
    interval: hb.ResourceFile | None = None,
    variant_types: list[Literal['INDEL', 'SNP', 'MNP', 'MIXED']] | None = None,
    job_attrs: dict | None = None,
    output_vcf_path: Path | None = None,
) -> Job:
    """
    Subset VCF to provided intervals.
    """
    if not interval and not variant_types:
        raise ValueError(
            'Either interval or variant_types must be defined for subset_vcf'
        )

    job_name = 'Subset VCF'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk SelectVariants'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    STANDARD.set_resources(j, ncpu=2)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    reference = fasta_res_group(b)
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    variant_types_param = ' '.join(
        f'--select-type-to-include {vt}' for vt in (variant_types or [])
    )
    cmd = f"""
    gatk SelectVariants \\
    -R {reference.base} \\
    -V {vcf['vcf.gz']} \\
    {f"-L {interval}" if interval else ''} \
    {variant_types_param} \\
    -O {j.output_vcf['vcf.gz']}
    """
    j.command(
        command(
            cmd,
            monitor_space=True,
        )
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: list[hb.ResourceFile],
    overwrite: bool = True,
    out_vcf_path: Path | None = None,
    site_only: bool = False,
    gvcf_count: int | None = None,
    job_attrs: dict | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`.

    Requires all VCFs to be strictly distinct, so doesn't work well
    for indels SelectVariants based on intervals from IntervalListTools,
    as ond indel might span 2 intervals and would end up in both.
    """
    if out_vcf_path and can_reuse(out_vcf_path, overwrite):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(out_vcf_path),
                'vcf.gz.tbi': f'{out_vcf_path}.tbi',
            }
        )

    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk GatherVcfsCloud'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))

    if gvcf_count:
        storage_gb = (1 if site_only else 2) * gvcf_count
        res = STANDARD.set_resources(j, fraction=1, storage_gb=storage_gb)
    else:
        res = STANDARD.set_resources(j, fraction=1)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v}' for v in input_vcfs])
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output $BATCH_TMPDIR/gathered.vcf.gz

    bcftools sort $BATCH_TMPDIR/gathered.vcf.gz -Oz \
    -o {j.output_vcf['vcf.gz']}
    
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """
    j.command(command(cmd, monitor_space=True))
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf
