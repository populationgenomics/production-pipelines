"""
Helper Hail Batch jobs useful for both individual and joint variant calling.
"""
from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group, command
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse


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
    out_vcf_path: Path | None = None,
    site_only: bool = False,
    sample_count: int | None = None,
    job_attrs: dict | None = None,
    sort: bool = False,
) -> tuple[list[Job], hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.

    Requires all VCFs to be strictly distinct, so doesn't work well
    for indels SelectVariants based on intervals from IntervalListTools,
    as ond indel might span 2 intervals and would end up in both.

    @param b: Batch object
    @param input_vcfs: list of Hail Batch ResourceFiles pointing to
        interval-split VCFs
    @param out_vcf_path: path to permanently write the resulting VCFs
    @param site_only: input VCFs are site-only
    @param sample_count: number of samples used for input VCFs (to determine the
        storage size)
    @param job_attrs: Hail Batch job attributes dictionary
    @param sort: whether to sort VCF before tabixing (computationally expensive, but
        required if the input VCFs can overlap)
    """
    if out_vcf_path and can_reuse([out_vcf_path, to_path(f'{out_vcf_path}.tbi')]):
        return [], b.read_input_group(
            **{
                'vcf.gz': str(out_vcf_path),
                'vcf.gz.tbi': f'{out_vcf_path}.tbi',
            }
        )

    jobs: list[Job] = []

    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk GatherVcfsCloud'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    if sample_count:
        storage_gb = (1 if site_only else 2) * sample_count
        res = STANDARD.request_resources(fraction=1, storage_gb=storage_gb)
    else:
        res = STANDARD.request_resources(fraction=1)
    res.set_to_job(j)

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
    --output {j.gathered_vcf}
    """
    j.command(command(cmd, monitor_space=True))
    jobs.append(j)

    job_attrs['tool'] = 'bcftools sort'
    j = b.new_job('Sort gathered VCF', job_attrs)
    j.image(image_path('bcftools'))
    if sample_count:
        storage_gb = (1 if site_only else 4) * sample_count
        STANDARD.set_resources(j, fraction=1, storage_gb=storage_gb)
    else:
        STANDARD.set_resources(j, fraction=1)
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = ''
    if sort:
        cmd += f"""
        bcftools sort {jobs[-1].gathered_vcf} \
        --temp-dir $BATCH_TMPDIR \
        -Oz -o {j.output_vcf['vcf.gz']}
        """

    cmd += f"""
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """

    j.command(command(cmd, monitor_space=True))
    jobs.append(j)
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))

    return jobs, j.output_vcf
