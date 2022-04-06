#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Loads the matrix table into ES.
"""

import os
import logging
import math
import subprocess
import click
import hail as hl
from hail_scripts.elasticsearch.hail_elasticsearch_client import HailElasticsearchClient

logger = logging.getLogger(__file__)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--es-host',
    'es_host',
    help='Elasticsearch host',
)
@click.option('--es-port', 'es_port', type=click.STRING, help='Elasticsearch port')
@click.option(
    '--es-username', 'es_username', type=click.STRING, help='Elasticsearch username'
)
@click.option(
    '--es-password', 'es_password', type=click.STRING, help='Elasticsearch password'
)
@click.option(
    '--es-index',
    'es_index',
    type=click.STRING,
    help='Elasticsearch index. Usually the dataset name. Will be lowercased',
    required=True,
)
@click.option(
    '--es-index-min-num-shards',
    'es_index_min_num_shards',
    default=1,
    help='Number of shards for the index will be the greater of this value '
    'and a calculated value based on the matrix.',
)
@click.option(
    '--use-spark', 'use_spark', is_flag=True, default=False,
)
def main(
    mt_path: str,
    es_host: str,
    es_port: str,
    es_username: str,
    es_password: str,
    es_index: str,
    es_index_min_num_shards: int,
    use_spark: bool,
):
    """
    Entry point.
    """
    if use_spark:
        hl.init(default_reference='GRCh38')
    else:
        from cpg_utils.hail import init_query_service
        init_query_service()

    if not all([es_host, es_port, es_username, es_password]):
        if any([es_host, es_port, es_username, es_password]):
            raise click.BadParameter(
                f'Either none, or all ES configuration parameters '
                f'must be specified: --es-host, --es-port, --es-username, --es-password. '
                f'If none are specified, defaults for the CPG are used'
            )
        es_host = 'elasticsearch.es.australia-southeast1.gcp.elastic-cloud.com'
        es_port = '9243'
        es_username = 'seqr'
        es_password = _read_es_password()

    es = HailElasticsearchClient(
        host=es_host,
        port=str(es_port),
        es_username=es_username,
        es_password=es_password,
        es_use_ssl=(es_host != 'localhost'),
    )

    mt = hl.read_matrix_table(mt_path)
    
    # Temporary fixes:
    mt = mt.annotate_rows(
        # AS_VQSLOD can be "Infinity" for indels , e.g.:
        # AS_VQSLOD=30.0692,18.2979,Infinity,17.5854,42.2131,1.5013
        # gs://cpg-seqr-main-tmp/seqr_loader/v0/AnnotateCohort/seqr_loader/checkpoints
        # /vqsr.ht
        # ht = hl.filter_intervals(ht, [hl.parse_locus_interval('chrX:52729395-52729396')])
        # hl.float() correctly parses this value, however, seqr loader doesn't 
        # recognise it, so we need to replace it with zero:
        info=mt.info.annotate(
            AS_VQSLOD=hl.if_else(
                hl.is_infinite(mt.info.AS_VQSLOD),
                0.0,
                mt.info.AS_VQSLOD,
            )
        ),
        # Seqr doesn't recognise PASS value in filters as pass, so need to remove it
        filters=mt.filters.filter(lambda val: val != 'PASS')
    )

    logger.info('Getting rows and exporting to the ES')
    row_table = elasticsearch_row(mt)
    es_shards = _mt_num_shards(mt, es_index_min_num_shards)

    es.export_table_to_elasticsearch(
        row_table,
        index_name=es_index.lower(),
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


def _mt_num_shards(mt, es_index_min_num_shards):
    """
    The greater of the user specified min shards and calculated based on the variants
    and samples
    """
    denominator = 1.4 * 10 ** 9
    calculated_num_shards = math.ceil((mt.count_rows() * mt.count_cols()) / denominator)
    return max(es_index_min_num_shards, calculated_num_shards)


def _cleanup(es, es_index, es_shards):
    # Current disk configuration requires the previous index to be deleted prior to large indices, ~1TB, transferring off loading nodes
    if es_shards < 25:
        es.wait_for_shard_transfer(es_index)


def _read_es_password(
    project_id='seqr-308602',
    secret_id='seqr-es-password',
    version_id='latest',
) -> str:
    """
    Read a GCP secret storing the ES password
    """
    password = os.environ.get('SEQR_ES_PASSWORD')
    if password:
        return password
    cmd = f'gcloud secrets versions access {version_id} --secret {secret_id} --project {project_id}'
    logger.info(cmd)
    return subprocess.check_output(cmd, shell=True).decode()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
