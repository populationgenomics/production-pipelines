#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils import to_path

from cpg_workflows.query_modules.seqr_loader import (
    subset_mt_to_samples,
    annotate_dataset_mt,
)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--sample-ids',
    'sample_ids_path',
    required=True,
)
@click.option(
    '--out-mt-path',
    'out_mt_path',
    required=True,
)
@click.option(
    '--checkpoint-prefix',
    'checkpoint_prefix',
    required=True,
)
def main(
    mt_path: str,
    sample_ids_path: str,
    out_mt_path: str,
    checkpoint_prefix: str,
):
    hl.init(default_reference='GRCh38')

    with to_path(sample_ids_path).open() as f:
        sample_ids = f.read().strip().split(',')

    subset_mt_path = to_path(checkpoint_prefix) / 'cohort-subset.mt'

    subset_mt_to_samples(
        mt_path=mt_path,
        sample_ids=sample_ids,
        out_mt_path=str(subset_mt_path),
    )

    annotate_dataset_mt(
        mt_path=str(subset_mt_path),
        out_mt_path=out_mt_path,
        checkpoint_prefix=checkpoint_prefix,
    )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
