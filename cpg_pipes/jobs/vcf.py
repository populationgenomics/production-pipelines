import logging
from typing import List, Tuple

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import resources, utils
from cpg_pipes.jobs import wrap_command

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: str = None,
    site_only: bool = False,
) -> Tuple[Job, hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    j = b.new_job(job_name)
    if utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j, b.read_input_group(**{
            'vcf.gz': output_vcf_path,
            'vcf.gz.tbi': output_vcf_path + '.tbi',
        })

    j.image(resources.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 3.75G/core ~ 7.5G
    j.storage(f'{1 + len(input_vcfs) * (0.1 if site_only else 2)}G')
        
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(wrap_command(f"""
    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{java_mem}g \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    """, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j, j.output_vcf
