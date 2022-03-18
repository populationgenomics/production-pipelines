"""
Helper Hail Batch jobs useful for both individual and joint variant calling.
"""

import logging
from typing import Tuple

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import images, buckets
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.storage import Path

logger = logging.getLogger(__file__)


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: list[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: Path | None = None,
    site_only: bool = False,
) -> Tuple[Job, hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    j = b.new_job(job_name)
    if output_vcf_path and buckets.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j, b.read_input_group(**{
            'vcf.gz': str(output_vcf_path),
            'vcf.gz.tbi': f'{output_vcf_path}.tbi',
        })

    j.image(images.GATK_IMAGE)
    STANDARD.set_resources(j, fraction=1)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(wrap_command(f"""
    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms25g \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    """, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf
