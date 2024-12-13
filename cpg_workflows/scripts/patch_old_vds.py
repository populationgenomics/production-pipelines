"""
Patch a VDS created prior to Hail 0.2.131 to include LGT calls in reference data.


If the sample is female (calculated by hail.impute_sex()), or the reference block is
within the PAR regions, the LGT is set to (0,0). Otherwise, it is set to (0)

Example usage:
    analysis-runner --dataset bioheart \
    --access-level standard \
    --output-dir vds-patch \
    --description 'Patch old VDS' \
    python3 patch_old_vds.py gs://cpg-bioheart-main/vds/tenk10k1-0.vds

"""

import argparse
from typing import TYPE_CHECKING

import hail as hl
from hail.vds.variant_dataset import VariantDataset

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
    tmp_output_path: str = output_path(f'{vds_name}-patch/reference_data', category='tmp')

    init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)
    vds: VariantDataset = hl.vds.read_vds(vds_path)

    # Patch the VDS to contain LGT calls as explained here: https://hail.zulipchat.com/#narrow/channel/123010-Hail-Query-0.2E2-support/topic/Error.20with.20multi_way_zip_join.20when.20combining.20two.20VDS
    chrx: VariantDataset = hl.vds.filter_chromosomes(vds, keep=['chrX'])
    biallelic_chrx: VariantDataset = chrx.variant_data.filter_rows(hl.len(chrx.variant_data.alleles) > 2, keep=False)
    imputed_sex: hl.Table = hl.impute_sex(biallelic_chrx.LGT, male_threshold=0.2)
    vds.reference_data = vds.reference_data.annotate_entries(
        LGT=hl.case()
        .when(~hl.is_defined(imputed_sex[vds.reference_data.s].is_female), hl.call(0, 0))
        .when(imputed_sex[vds.reference_data.s].is_female, hl.call(0, 0))
        .when(vds.reference_data.locus.in_autosome_or_par(), hl.call(0, 0))
        .default(hl.call(0)),
    )
    vds.reference_data.write(tmp_output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vds-path', help='Path to VDS to patch reference data.')
    args: argparse.Namespace = parser.parse_args()
    main(args.vds_path)
