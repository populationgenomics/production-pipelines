#!/usr/bin/env python3

"""
Run VEP in parallel using batch backend
"""

import click

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path, output_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.jobs.vep import add_vep_jobs
from cpg_workflows.query_modules import seqr_loader


@click.command()
@click.option('--vcf-path', required=True, help='Full path to VCF to run VEP on')
@click.option(
    '--output-ht',
    required=True,
    help='Path to where finished VEP-annotated VCF will be output',
)
@click.option(
    '--to-mt',
    is_flag=True,
    help='Complete transition to MT, only available on HT output',
)
def main(vcf_path: str, output_ht: str, to_mt: bool = False):
    """
    Run VEP in parallel using Picard tools intervals as partitions.
    Input: the full path to a VCF file, along with a tabix (.tbi) file,
    located in the same directory.
    """
    vep_version = get_config()['workflow']['vep_version']
    vep_image = image_path(f'vep_{vep_version}')
    scatter_count = get_config()['vep']['scatter_count']
    b = get_batch(f'Run VEP with Batch Backend, image {vep_image}, scatter count {scatter_count}')
    vep_ht = output_path(output_ht)
    vep_jobs = add_vep_jobs(
        b=b,
        input_vcfs=[to_path(vcf_path)],
        tmp_prefix=to_path(output_path('vcf_fragments/', 'tmp')),
        out_path=to_path(vep_ht),
        scatter_count=scatter_count,
    )
    if to_mt:
        assert vep_ht.endswith('.ht')
        vep_mt = vep_ht[: -len('.ht')] + '.mt'
        j = b.new_job('annotate cohort', {'tool': 'hail query'})
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                seqr_loader,
                seqr_loader.annotate_cohort.__name__,
                str(vcf_path),
                vep_mt,
                vep_ht,
                setup_gcp=True,
            ),
        )
        if vep_jobs:
            j.depends_on(*vep_jobs)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
