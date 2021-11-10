#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster. 

Loads the matrix table into ES.
"""

import logging
import math
import subprocess
import click
import hail as hl

from hail_scripts.v02.utils.elasticsearch_client import ElasticsearchClient

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--mt-path',
    'mt_path',
    required=True,
)
@click.option(
    '--bucket',
    'work_bucket',
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
    '--prod', is_flag=True, help='Run under the production ES credentials instead'
)
def main(
    mt_path: str,
    es_host: str,
    es_port: str,
    es_username: str,
    es_password: str,
    es_index: str,
    es_index_min_num_shards: int,
    prod: bool,  # pylint: disable=unused-argument
):  # pylint: disable=missing-function-docstring
    hl.init(default_reference='GRCh38')

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

    es = ElasticsearchClient(
        host=es_host,
        port=str(es_port),
        es_username=es_username,
        es_password=es_password,
        es_use_ssl=(es_host != 'localhost'),
    )

    mt = hl.read_matrix_table(mt_path)

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
    """
    # Converts a mt to the row equivalent.
    ht = mt.rows()
    ht.describe()
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    table = ht.drop('vep').flatten()
    # When flattening, the table is unkeyed, which causes problems because our locus and alleles should not
    # be normal fields. We can also re-key, but I believe this is computational?
    table = table.drop(table.locus, table.alleles)
    return table


def _read_es_password(
    project_id='seqr-308602',
    secret_id='seqr-es-password',
    version_id='latest',
) -> str:
    """
    Read a GCP secret storing the ES password
    """
    cmd = f'gcloud secrets versions access {version_id} --secret {secret_id} --project {project_id}'
    logger.info(cmd)
    return subprocess.check_output(cmd, shell=True).decode()


def _mt_num_shards(mt, es_index_min_num_shards):
    """
    The greater of the user specified min shards and calculated based on the variants
    and samples
    """
    denominator = 1.4 * 10 ** 9
    calculated_num_shards = math.ceil((mt.count_rows() * mt.count_cols()) / denominator)
    return max(es_index_min_num_shards, calculated_num_shards)


def _cleanup(es, es_index, es_shards):
    es.route_index_off_temp_es_cluster(es_index)
    # Current disk configuration requires the previous index to be deleted prior to large indices, ~1TB, transferring off loading nodes
    if es_shards < 25:
        es.wait_for_shard_transfer(es_index)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
