"""Patch a VDS created prior to Hail 0.2.131 to include LGT calls in reference data."""

import argparse
from typing import TYPE_CHECKING

import hail as hl

from cpg_utils.config import output_path
from cpg_utils.hail_batch import init_batch

if TYPE_CHECKING:
    from hail.vds.variant_dataset import VariantDataset


def main(vds_path: str) -> None:
    """Patch the VDS and write to a tmp bucket.

    Args:
        vds_path (str): Full gs:// path to the VDS to patch
    """
    vds_name: str = vds_path.split('/')[-1].replace('.vds', '')
    tmp_output_path: str = output_path(vds_name, category='tmp')

    init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)
    vds: VariantDataset = hl.vds.read_vds(vds_path)

    # Patch the VDS to contain LGT calls as explained here: https://hail.zulipchat.com/#narrow/channel/123010-Hail-Query-0.2E2-support/topic/Error.20with.20multi_way_zip_join.20when.20combining.20two.20VDS
    vds.reference_data = vds.reference_data.annotate_entries(
        LGT=hl.if_else(hl.Call.is_haploid, hl.call(0), hl.call(0, 0)),
    )
    vds.reference_data.write(tmp_output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vds-path', help='Path to VDS to patch reference data.')
    args: argparse.Namespace = parser.parse_args()
    main(args.vds_path)
