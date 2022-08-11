"""
Helper Hail Batch jobs useful for both individual and joint variant calling.
"""

import logging

import hailtop.batch as hb
from cpg_utils.hail_batch import image_path, fasta_res_group
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD

logger = logging.getLogger(__file__)


def subset_vcf(
    b: hb.Batch,
    vcf: hb.ResourceGroup,
    interval: hb.Resource,
    job_attrs: dict | None = None,
    output_vcf_path: Path | None = None,
) -> Job:
    """
    Subset VCF to provided intervals.
    """
    job_name = 'Subset VCF'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    STANDARD.set_resources(j, ncpu=2)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    reference = fasta_res_group(b)

    cmd = f"""
    gatk SelectVariants \\
    -R {reference.base} \\
    -V {vcf['vcf.gz']} \\
    -L {interval} \\
    -O {j.output_vcf['vcf.gz']}
    """
    j.command(
        wrap_command(
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
) -> tuple[list[Job], hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    if out_vcf_path and utils.can_reuse(out_vcf_path, overwrite):
        return [], b.read_input_group(
            **{
                'vcf.gz': str(out_vcf_path),
                'vcf.gz.tbi': f'{out_vcf_path}.tbi',
            }
        )

    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
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
    --output {j.output_vcf['vcf.gz']}
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))
    return [j], j.output_vcf
