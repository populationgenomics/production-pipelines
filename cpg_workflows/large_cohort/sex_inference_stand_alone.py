#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs annotate_sex()

 analysis-runner --dataset "bioheart" \
    --description "standalone annotate-sex method" \
    --access-level "test" \
    --output-dir "qc-stand-alone/annotate-sex" \
    sex_inference_stand_alone.py --input-dir=gs://cpg-bioheart-test/vds/5-0.vds

"""

import logging

import hail as hl
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, genome_build, output_path, init_batch
from gnomad.sample_qc.pipeline import annotate_sex
from hail.vds.variant_dataset import VariantDataset

import click


def impute_sex(
    vds_path: str,
) -> hl.Table:
    """
    Impute sex based on coverage.
    """

    init_batch()

    # run()
    vds = hl.vds.read_vds(str(vds_path))

    # Remove centromeres and telomeres:
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    if tel_cent_ht.count() > 0:
        vds = hl.vds.filter_intervals(vds, tel_cent_ht, keep=False)

    # Load calling intervals
    calling_intervals_path = reference_path(f'broad/genome_calling_interval_lists')
    calling_intervals_ht = hl.import_locus_intervals(
        str(calling_intervals_path), reference_genome=genome_build()
    )
    logging.info('Calling intervals table:')
    calling_intervals_ht.describe()

    # Pre-filter here and setting `variants_filter_lcr` and `variants_filter_segdup`
    # below to `False` to avoid the function calling gnomAD's `resources` module:
    for name in ['lcr_intervals_ht', 'seg_dup_intervals_ht']:
        interval_table = hl.read_table(str(reference_path(f'gnomad/{name}')))
        if interval_table.count() > 0:
            # remove all rows where the locus falls within a defined interval
            tmp_variant_data = vds.variant_data.filter_rows(
                hl.is_defined(interval_table[vds.variant_data.locus]), keep=False
            )
            vds = VariantDataset(
                reference_data=vds.reference_data, variant_data=tmp_variant_data
            ).checkpoint(str(f'{name}_checkpoint.vds'))
            logging.info(f'count post {name} filter:{vds.variant_data.count()}')

    # Infer sex (adds row fields: is_female, var_data_chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        tmp_prefix=str('annotate_sex'),
        overwrite=not get_config()['workflow'].get('check_intermediates'),
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
        variants_only_x_ploidy=True,
        variants_only_y_ploidy=False,
        variants_filter_lcr=False,  # already filtered above
        variants_filter_segdup=False,  # already filtered above
        variants_filter_decoy=False,
    )
    logging.info('Sex table:')
    sex_ht.describe()
    sex_ht = sex_ht.transmute(
        impute_sex_stats=hl.struct(
            f_stat=sex_ht.f_stat,
            n_called=sex_ht.n_called,
            expected_homs=sex_ht.expected_homs,
            observed_homs=sex_ht.observed_homs,
        )
    )

    # output writing
    sex_ht.write(output_path(f'sex.ht', 'analysis'))


@click.option(
    '--vds-path',
    help='GCS file path to VDS',
    type=str,
)
@click.command()
def main(vds_path):
    impute_sex(vds_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
