"""
Batch script that runs a VEP job on input VCFs.
"""
import os.path
from typing import List

from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.jobs import vep
import click


@click.command()
@click.option(
    '--vcf',
    '--vcf-path',
    'vcf_paths',
    multiple=True,
    help='Paths to VCFs (can be multiple)',
)
@click.option(
    '--project',
    'project',
    help='Billing project for Batch'
)
def main(vcf_paths: List[str], project: str):
    b = setup_batch('Run VEP', billing_project=project)

    for vcf_path in vcf_paths:
        out_vcf_path = os.path.splitext(vcf_path)[0] + '-vep.vcf.bgz'
        vep.vep(
            b, 
            vcf_path=vcf_path, 
            out_vcf_path=out_vcf_path,
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
