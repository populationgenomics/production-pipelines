#!/usr/bin/env python3

"""
Dataproc script to annotate cohort.
"""

import click
import hail as hl

from cpg_utils import to_path

from cpg_workflows.query_modules.seqr_loader import subset_mt_to_samples
from cpg_workflows.query_modules.seqr_loader_sv import annotate_dataset_sv


@click.command()
@click.option('--mt-path', required=True)
@click.option('--sgids', required=True)
@click.option('--out-mt-path', required=True)
@click.option('--checkpoint-prefix', required=True)
def main(mt_path: str, sgids: str, out_mt_path: str, checkpoint_prefix: str):
    """
    Schedule two job stages - subset_mt_to_samples and annotate_dataset_sv
    Args:
        mt_path (str): path to the reformatted-annotations MT
        sgids (str): path to a file containing subset of sample IDs
        out_mt_path (str): where to write the final MT
        checkpoint_prefix ():
    """

    hl.init(default_reference='GRCh38')

    with to_path(sgids).open() as f:
        sgid_list = f.read().strip().split(',')

    subset_mt = str(to_path(checkpoint_prefix) / 'cohort-subset.mt')

    subset_mt_to_samples(mt_path=mt_path, sample_ids=sgid_list, out_mt_path=subset_mt)

    annotate_dataset_sv(mt_path=subset_mt, out_mt_path=out_mt_path)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
