"""
Helper Hail Batch jobs useful for both individual and joint variant calling.
"""
from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group, command
from cpg_workflows.resources import STANDARD, storage_for_joint_vcf
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
) -> list[Job]:
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
    @param sort: whether to sort VCF before tabixing (computationally expensive,
        but required if the input VCFs can overlap)
    """
    jobs: list[Job | None] = []
    gathered_vcf: hb.ResourceFile

    if not can_reuse(out_vcf_path):
        job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
        j = b.new_job(job_name, (job_attrs or {}) | {'tool': 'gatk GatherVcfsCloud'})
        j.image(image_path('gatk'))
        res = STANDARD.set_resources(
            j, storage_gb=storage_for_joint_vcf(sample_count, site_only)
        )

        input_cmdl = ' '.join([f'--input {v}' for v in input_vcfs])
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
        mv $BATCH_TMPDIR/gathered.vcf.gz {j.output_vcf}
        """
        j.command(command(cmd, monitor_space=True))
        if not sort and out_vcf_path:
            b.write_output(j.output_vcf, str(out_vcf_path))
        assert isinstance(j.output_vcf, hb.ResourceFile)
        gathered_vcf = j.output_vcf
    else:
        gathered_vcf = b.read_input(str(out_vcf_path))

    out_tbi_path = to_path(f'{out_vcf_path}.tbi')
    if not can_reuse(out_tbi_path):
        if sort:
            jobs.append(
                sort_vcf(
                    b,
                    vcf=gathered_vcf,
                    out_vcf_path=out_vcf_path,
                    job_attrs=job_attrs,
                    sample_count=sample_count,
                    site_only=site_only,
                )
            )
        else:
            jobs.append(
                tabix_vcf(
                    b,
                    vcf=gathered_vcf,
                    out_tbi_path=out_tbi_path,
                    job_attrs=job_attrs,
                    sample_count=sample_count,
                    site_only=site_only,
                )
            )
    return [j for j in jobs if j is not None]


def sort_vcf(
    b: hb.Batch,
    vcf: hb.ResourceFile,
    out_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
    site_only: bool = False,
    sample_count: int | None = None,
) -> Job | None:
    """
    Sort and index VCF.
    """
    if can_reuse([out_vcf_path, f'{out_vcf_path}.tbi']):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'bcftools sort'}
    j = b.new_job('Sort gathered VCF', job_attrs)
    j.image(image_path('bcftools'))
    if storage_gb := storage_for_joint_vcf(sample_count, site_only):
        storage_gb *= 2  # sort needs extra tmp space
    STANDARD.set_resources(j, fraction=1, storage_gb=storage_gb)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    bcftools sort {vcf} \
    --temp-dir $BATCH_TMPDIR \
    -Oz -o {j.output_vcf['vcf.gz']}

    tabix -p vcf {j.output_vcf['vcf.gz']}
    """

    j.command(command(cmd, monitor_space=True))
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))

    return j


def tabix_vcf(
    b: hb.Batch,
    vcf: hb.ResourceFile,
    out_tbi_path: Path | None = None,
    job_attrs: dict | None = None,
    site_only: bool = False,
    sample_count: int | None = None,
) -> Job | None:
    """
    Index VCF, assuming it's sorted.
    """
    if can_reuse(out_tbi_path):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'tabix'}
    j = b.new_job('Tabix sorted VCF', job_attrs)
    j.image(image_path('bcftools'))
    STANDARD.set_resources(
        j, fraction=1, storage_gb=storage_for_joint_vcf(sample_count, site_only)
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    mv {vcf} $BATCH_TMPDIR/result.vcf.gz
    tabix -p vcf $BATCH_TMPDIR/result.vcf.gz
    mv $BATCH_TMPDIR/result.vcf.gz.tbi {j.output_tbi}
    """

    j.command(command(cmd, monitor_space=True))
    if out_tbi_path:
        b.write_output(j.output_tbi, str(out_tbi_path))

    return j
