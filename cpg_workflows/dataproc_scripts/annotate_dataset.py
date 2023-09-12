#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils import to_path

from cpg_workflows.query_modules.seqr_loader import (
    subset_mt_to_sgids,
    annotate_dataset_mt,
)


@click.command()
@click.option('--mt-path', required=True)
@click.option('--sgids', required=True)
@click.option('--out-mt-path', required=True)
@click.option('--checkpoint-prefix', required=True)
def main(mt_path: str, sgids: str, out_mt_path: str, checkpoint_prefix: str):
    hl.init(default_reference='GRCh38')

    with to_path(sgids).open() as f:
        sgid_list = f.read().strip().split(',')

    subset_mt_path = to_path(checkpoint_prefix) / 'cohort-subset.mt'

    subset_mt_to_sgids(
        mt_path=mt_path, sgid_list=sgid_list, out_mt_path=str(subset_mt_path)
    )

    annotate_dataset_mt(
        mt_path=str(subset_mt_path),
        out_mt_path=out_mt_path,
        checkpoint_prefix=checkpoint_prefix,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
