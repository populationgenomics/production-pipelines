#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Loads the matrix table into an ElasticSearch index.
"""

import logging
import math
import click
import hail as hl
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config
from hail_scripts.elasticsearch.hail_elasticsearch_client import HailElasticsearchClient

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--es-index',
    'es_index',
    type=click.STRING,
    help='Elasticsearch index. Usually the dataset name. Will be lowercased',
    required=True,
)
@click.option(
    '--use-spark',
    'use_spark',
    is_flag=True,
    default=False,
)
@click.option(
    '--liftover-path',
    'liftover_path',
    help='Path to liftover chain',
    required=True,
)
@click.option(
    '--es-password',
    'password',
)
def main(
    mt_path: str,
    es_index: str,
    use_spark: bool,
    liftover_path: str,
    password: str = None,
):
    """
    Entry point.
    """
    es_index = es_index.lower()
    
    if use_spark:
        hl.init(default_reference='GRCh38')
    else:
        from cpg_utils.hail_batch import init_batch
        init_batch()
        
    host = get_config()['elasticsearch']['host']
    port = str(get_config()['elasticsearch']['port'])
    username = get_config()['elasticsearch']['username']
    password = password or read_secret(
        project_id=get_config()['elasticsearch']['password_project_id'],
        secret_name=get_config()['elasticsearch']['password_secret_id'],
        fail_gracefully=False,
    )
    assert password

    es = HailElasticsearchClient(
        host=host,
        port=port,
        es_username=username,
        es_password=password,
        es_use_ssl=(host != 'localhost'),
    )

    mt = hl.read_matrix_table(mt_path)
    # Annotate GRCh37 coordinates here, as they are not supported by Batch Backend
    logger.info('Adding GRCh37 coords')
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_path, rg37)
    mt = mt.annotate_rows(rg37_locus=hl.liftover(mt.locus, 'GRCh37'))

    logger.info('Getting rows and exporting to the ES')
    row_table = elasticsearch_row(mt)
    es_shards = _mt_num_shards(mt)

    es.export_table_to_elasticsearch(
        row_table,
        index_name=es_index,
        num_shards=es_shards,
        write_null_values=True,
    )
    _cleanup(es, es_index, es_shards)


def elasticsearch_row(mt: hl.MatrixTable):
    """
    Prepares the mt to export using ElasticsearchClient V02.
    - Flattens nested structs
    - drops locus and alleles key

    Borrowed from:
    https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/495f0d1b4d49542557ca5cccf98a23fc627260bf/luigi_pipeline/lib/model/seqr_mt_schema.py
    """
    # Converts a mt to the row equivalent.
    ht = mt.rows()
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    table = ht.drop('vep').flatten()
    # When flattening, the table is unkeyed, which causes problems because our locus
    # and alleles should not
    # be normal fields. We can also re-key, but I believe this is computational?
    table = table.drop(table.locus, table.alleles)
    table.describe()
    return table


def _mt_num_shards(mt):
    """
    Calculate number of shareds from the number of variants and samples.
    """
    denominator = 1.4 * 10**9
    calculated_num_shards = math.ceil((mt.count_rows() * mt.count_cols()) / denominator)
    return calculated_num_shards


def _cleanup(es, es_index, es_shards):
    # Current disk configuration requires the previous index to be deleted prior to large indices, ~1TB, transferring off loading nodes
    if es_shards < 25:
        es.wait_for_shard_transfer(es_index)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
