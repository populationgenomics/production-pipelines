#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils.hail_batch import genome_build
from cpg_workflows.query_modules.seqr_loader import annotate_cohort


@click.command()
@click.option(
    '--vcf-path',
    'vcf_path',
    required=True,
)
@click.option(
    '--out-mt-path',
    'out_mt_path',
    required=True,
)
@click.option(
    '--vep-ht-path',
    'vep_ht_path',
    required=True,
)
@click.option(
    '--siteonly-vqsr-vcf-path',
    'siteonly_vqsr_vcf_path',
)
@click.option(
    '--checkpoint-prefix',
    'checkpoint_prefix',
)
def main(
    vcf_path: str,
    out_mt_path: str,
    vep_ht_path: str,
    siteonly_vqsr_vcf_path: str | None = None,
    checkpoint_prefix: str | None = None,
):
    hl.init(default_reference=genome_build())

    annotate_cohort(
        vcf_path=vcf_path,
        out_mt_path=out_mt_path,
        vep_ht_path=vep_ht_path,
        site_only_vqsr_vcf_path=siteonly_vqsr_vcf_path,
        checkpoint_prefix=checkpoint_prefix,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
