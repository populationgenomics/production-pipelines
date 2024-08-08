"""
This Stage does some Hail Query work to add background and Metadata to the QC tables
"""

import logging
from argparse import ArgumentParser
from os.path import join

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch


def prepare_background_mt(data_path: str, dataset: str, qc_variants_ht: hl.Table, tmp_path) -> hl.MatrixTable:
    """

    Args:
        data_path ():
        qc_variants_ht ():
        tmp_path ():

    Returns:
        A densified MT, checkpointed
    """
    background_mt_checkpoint_path = join(tmp_path, f'{dataset}_densified_background.mt')

    if to_path(data_path).suffix == '.mt':
        background_mt = hl.read_matrix_table(data_path)
        background_mt = hl.split_multi(background_mt, filter_changed_loci=True)
        background_mt = background_mt.semi_join_rows(qc_variants_ht)
        return background_mt.densify().checkpoint(background_mt_checkpoint_path, overwrite=True)
    elif to_path(data_path).suffix == '.vds':
        background_vds = hl.vds.read_vds(data_path)
        background_vds = hl.vds.split_multi(background_vds, filter_changed_loci=True)
        background_vds = hl.vds.filter_variants(background_vds, qc_variants_ht)
        background_mt = hl.vds.to_dense_mt(background_vds)
        logging.info(f'Checkpointing background_mt to {background_mt_checkpoint_path}')
        return background_mt.checkpoint(background_mt_checkpoint_path, overwrite=True)
    else:
        raise ValueError('Background dataset path must be either .mt or .vds')


def prepare_metadata_table(background_sample_qc_paths: list[str], sample_qc_ht: hl.Table) -> hl.Table:
    """

    Args:
        background_sample_qc_paths ():
        sample_qc_ht (hl.Table):

    Returns:
        the union of all Metadata tables, with columns reordered WRT the SampleQC Hail Table
    """
    # annotate background mt with metadata info derived from SampleQC stage
    allow_missing_columns = config_retrieve(['large_cohort', 'pca_background', 'allow_missing_columns'], False)
    metadata_tables = []
    for path in background_sample_qc_paths:
        metadata_tables.append(hl.read_table(path))

    metadata_table = hl.Table.union(*metadata_tables, unify=allow_missing_columns)

    # reorder columns
    reference_fields = sample_qc_ht.sample_qc.dtype.fields
    ordered_sample_qc = {field: metadata_table.sample_qc[field] for field in reference_fields}
    return metadata_table.annotate(sample_qc=hl.struct(**ordered_sample_qc))


def add_background(dense_mt: hl.MatrixTable, sample_qc_ht: hl.Table, tmp_path: str) -> tuple[hl.MatrixTable, hl.Table]:
    """
    Add background dataset samples to the dense MT and sample QC HT.

    Args:
        dense_mt ():
        sample_qc_ht ():
        tmp_path (): where to generate temp data, if required

    Returns:
        Annotated versions of the Sample QC HT, and Dense MT
    """
    pca_bg_config = config_retrieve(['large_cohort', 'pca_background'])

    # whole table read in from config path
    qc_variants_ht = hl.read_table(config_retrieve(['references', 'ancestry', 'sites_table']))

    allow_missing_columns = pca_bg_config.get('allow_missing_columns', False)

    # loop over each dataset. If this is empty, we don't run this stage
    # happy to hard pull from config and fail if absent
    for dataset in pca_bg_config['datasets']:

        dataset_dict = pca_bg_config[dataset]
        logging.info(f'Adding background dataset {dataset_dict["dataset_path"]}')

        # prepare the background MT, from either a Mt or VDS
        background_mt = prepare_background_mt(dataset_dict['dataset_path'], dataset, qc_variants_ht, tmp_path=tmp_path)
        metadata_table = prepare_metadata_table(dataset_dict['metadata_table'], sample_qc_ht)

        background_mt = background_mt.annotate_cols(**metadata_table[background_mt.col_key])

        if populations_to_filter := pca_bg_config.get('superpopulation_to_filter', False):
            logging.info(f'Filtering background samples by {populations_to_filter}')
            background_mt = background_mt.filter_cols(
                hl.literal(populations_to_filter).contains(background_mt.superpopulation),
            )
            logging.info(f'Finished filtering background, kept samples that are {populations_to_filter}')
        if background_relateds_to_drop := config_retrieve(
            ['large_cohort', 'pca_background', dataset, 'background_relateds_to_drop'],
            False,
        ):
            logging.info(
                f'Removing related samples from background dataset {dataset}.\n'
                f'Background relateds to drop: {background_relateds_to_drop}',
            )
            background_relateds_to_drop_ht = hl.read_table(background_relateds_to_drop)
            background_mt = background_mt.filter_cols(
                ~hl.is_defined(background_relateds_to_drop_ht[background_mt.col_key]),
            )
        else:
            logging.info('No related samples to drop from background dataset')

        # save metadata info before merging dense and background datasets
        sample_qc_ht = sample_qc_ht.union(background_mt.cols(), unify=allow_missing_columns)

        background_mt = (
            background_mt.select_cols().select_rows().select_entries('GT', 'GQ', 'DP', 'AD').naive_coalesce(5000)
        )
        # combine dense dataset with background population dataset
        dense_mt = dense_mt.union_cols(background_mt)

    if drop_columns := pca_bg_config.get('drop_columns'):
        sample_qc_ht = sample_qc_ht.drop(*drop_columns)

    return dense_mt, sample_qc_ht


def cli_main():
    """
    A command-line entrypoint for the ancestry add-background process
    """

    parser = ArgumentParser()
    parser.add_argument('--qc_in', help='The QC Table to read in')
    parser.add_argument('--dense_in', help='The Dense MT to read in')
    parser.add_argument('--qc_out', help='The output path for the QC table')
    parser.add_argument('--dense_out', help='The output path for the dense MT')
    parser.add_argument('--tmp', help='Path to write temporary data to')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    main(qc_in=args.qc_in, dense_in=args.dense_in, qc_out=args.qc_out, dense_out=args.dense_out, tmp_path=args.tmp)


def main(qc_in: str, dense_in: str, qc_out: str, dense_out: str, tmp_path: str):
    """

    Args:
        qc_in (str):
        dense_in (str):
        qc_out (str):
        dense_out (str):
        tmp_path (str):
    """

    # start a Hail Query (on-batch) Runtime
    init_batch()

    dense_mt = hl.read_matrix_table(dense_in).select_entries('GT', 'GQ', 'DP', 'AD')
    sample_qc_ht = hl.read_table(qc_in)

    dense_mt, sample_qc_ht = add_background(dense_mt, sample_qc_ht, tmp_path)

    logging.info(f'Writing dense_mt to {dense_out}')
    dense_mt.write(dense_out)

    logging.info(f'Writing sample_qc_ht to {qc_out}')
    sample_qc_ht.write(qc_out, overwrite=True)


if __name__ == '__main__':
    cli_main()
