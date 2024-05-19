#!/usr/bin/env python3

"""
Hail script to submit on a dataproc cluster.

Loads the matrix table into an ElasticSearch index.
"""

import logging
import math
from pprint import pformat

import click
import coloredlogs
import elasticsearch

import hail as hl

from cpg_utils import to_path
from cpg_utils.cloud import read_secret
from cpg_utils.config import get_config, reference_path
from hail_scripts.elasticsearch.hail_elasticsearch_client import HailElasticsearchClient

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='INFO', fmt=fmt)


class HailElasticsearchClientV8(HailElasticsearchClient):
    """
    The Broad's seqr-loading-pipelines pins the Elasticsearch client to v7. We use v8,
    so we are overriding the class to adjust for Elasticsearch v8.
    """

    # noinspection PyMissingConstructor
    def __init__(
        self,
        host: str,
        port: str,
        es_username: str,
        es_password: str,
    ):
        """
        Overriding base __init__: the difference with v7 is that in v8,
        elasticsearch.Elasticsearch constructor takes one URL string, in contrast
        with v7 which takes "host" and "port" strings separately. Note that we are
        not calling the base init, which would fail on v8. Instead, we are completely
        overriding init.
        """
        self._host = host
        self._port = port
        self._es_username = es_username
        self._es_password = es_password
        self._es_use_ssl = True

        auth = (self._es_username, self._es_password) if self._es_password else None

        if not host.startswith('http://') or not host.startswith('https://'):
            scheme = 'https' if self._es_use_ssl else 'http'
            host = f'{scheme}://{host}'
        _host = f'{host}:{port}'
        self.es = elasticsearch.Elasticsearch(_host, basic_auth=auth)

        # check connection
        logging.info(pformat(self.es.info()))

    def export_table_to_elasticsearch(self, *args, **kwargs):
        """Override to adjust for ES V7."""

        # Copied from older hail_scripts/v02/utils/elasticsearch_client.py:
        # https://github.com/populationgenomics/hail-elasticsearch-pipelines/blob/main/hail_scripts/v02/utils/elasticsearch_client.py#L128-L133
        # that's before elasticsearch in upstream seqr-loading-pipelines was pinned to v7.
        # Without this config, ES API errors with the following:
        # > Cannot detect ES version - typically this happens if the network/Elasticsearch
        # cluster is not accessible or when targeting a WAN/Cloud instance without the
        # proper setting 'es.nodes.wan.only'
        if self._es_use_ssl:
            kwargs.setdefault('elasticsearch_config', {}).update(
                {
                    'es.net.ssl': 'true',
                    # If using SSL, the instance is likely managed, in which case we
                    # can't discover nodes.
                    'es.nodes.wan.only': 'true',
                },
            )

        # deprecated in ES=V8:
        kwargs['index_type_name'] = ''

        super().export_table_to_elasticsearch(*args, **kwargs)


@click.command()
@click.option('--mt-path', required=True)
@click.option(
    '--es-index',
    help='Elasticsearch index. Usually the dataset name. Will be lowercased',
    required=True,
)
@click.option('--done-flag-path', help='File to touch in the end')
@click.option('--es-password')
def main(
    mt_path: str,
    es_index: str,
    done_flag_path: str,
    es_password: str | None = None,
):
    """
    Entry point.
    """
    es_index = es_index.lower()
    hl.init(default_reference='GRCh38')

    host = get_config()['elasticsearch']['host']
    port = str(get_config()['elasticsearch']['port'])
    username = get_config()['elasticsearch']['username']
    project_id = get_config()['elasticsearch']['password_project_id']
    secret_name = get_config()['elasticsearch']['password_secret_id']
    password = es_password or read_secret(
        project_id=project_id,
        secret_name=secret_name,
        fail_gracefully=False,
    )
    assert password

    print(f'Connecting to ElasticSearch: host="{host}", port="{port}", user="{username}"')
    print(f'Reading passport from secret "{secret_name}" in project "{project_id}"')
    es = HailElasticsearchClientV8(
        host=host,
        port=port,
        es_username=username,
        es_password=password,
    )

    mt = hl.read_matrix_table(mt_path)

    logging.info('Getting rows and exporting to the ES')
    row_ht = elasticsearch_row(mt)
    es_shards = _mt_num_shards(mt)

    es.export_table_to_elasticsearch(
        row_ht,
        index_name=es_index,
        num_shards=es_shards,
        write_null_values=True,
    )
    _cleanup(es, es_index, es_shards)
    with to_path(done_flag_path).open('w') as f:
        f.write('done')


def elasticsearch_row(mt: hl.MatrixTable):
    """
    Prepares the mt for export.
    - Flattens nested structs
    - drops locus and alleles key
    Borrowed from:
    https://github.com/broadinstitute/hail-elasticsearch-pipelines/blob/main/luigi_pipeline/lib/model/seqr_mt_schema.py
    """
    # Converts a mt to the row equivalent.
    if isinstance(mt, hl.MatrixTable):
        mt = mt.rows()
    if 'vep' in mt.row:
        mt = mt.drop('vep')
    key = mt.key
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    table = mt.flatten()
    # When flattening, the table is unkeyed, which causes problems because our row keys should not
    # be normal fields. We can also re-key, but I believe this is computational?
    # PS: row key is often locus and allele, but does not have to be
    table = table.drop(*key)
    table.describe()
    return table


def _mt_num_shards(mt):
    """
    Calculate the number of shards from the number of variants and sequencing groups.
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
