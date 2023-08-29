#!/usr/bin/env python3

"""
Run VEP in parallel using batch backend
"""

import click

from cpg_utils import to_path
from cpg_utils.hail_batch import dataset_path, output_path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_workflows.jobs.vep import add_vep_jobs


@click.command()
@click.option('--vcf-path', required=True, help='Full path to VCF to run VEP on')
@click.option(
    '--output-ht',
    required=True,
    help='Path to where finished VEP-annotated VCF will be output',
)
def main(vcf_path: str, output_ht: str):
    """
    Run VEP in parallel using Picard tools intervals as partitions.
    Input: the full path to a VCF file, along with a tabix (.tbi) file,
    located in the same directory.
    """
    vep_image = get_config()['images']['vep']
    b = get_batch(f'Run VEP with Batch Backend, image {vep_image}')
    add_vep_jobs(
        b=b,
        input_siteonly_vcf_path=to_path(vcf_path),
        tmp_prefix=to_path(output_path('vcf_fragments/', 'tmp')),
        out_path=to_path(dataset_path(output_ht)),
        scatter_count=get_config()['vep']['scatter_count'],
    )
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
