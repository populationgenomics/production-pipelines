#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils.hail_batch import genome_build
from cpg_workflows.query_modules.seqr_loader_sv import annotate_cohort_sv


@click.command()
@click.option('--vcf-path', 'vcf_path', required=True)
@click.option('--out-mt-path', 'out_mt_path', required=True)
@click.option('--checkpoint-prefix', 'checkpoint_prefix')
def main(vcf_path: str, out_mt_path: str, checkpoint_prefix: str | None = None):
    hl.init(default_reference=genome_build())

    annotate_cohort_sv(
        vcf_path=vcf_path,
        out_mt_path=out_mt_path,
        checkpoint_prefix=checkpoint_prefix,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
