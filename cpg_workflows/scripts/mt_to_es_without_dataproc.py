"""
I've peppered this with references to the original source code I'm borrowing/stealing from

major sources:
https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/hail_elasticsearch_client.py
https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/luigi_pipeline/lib/hail_tasks.py
"""

import logging
import math
import time
from argparse import ArgumentParser
from io import StringIO
from sys import exit

import elasticsearch

import hail as hl

from cpg_utils import to_path
from cpg_utils.cloud import read_secret
from cpg_utils.config import config_retrieve

# CONSTANTS stolen from https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_utils.py#L13
# make encoded values as human-readable as possible
ES_FIELD_NAME_ESCAPE_CHAR = '$'
ES_FIELD_NAME_BAD_LEADING_CHARS = {'_', '-', '+', ES_FIELD_NAME_ESCAPE_CHAR}
ES_FIELD_NAME_SPECIAL_CHAR_MAP = {
    '.': '_$dot$_',
    ',': '_$comma$_',
    '#': '_$hash$_',
    '*': '_$star$_',
    '(': '_$lp$_',
    ')': '_$rp$_',
    '[': '_$lsb$_',
    ']': '_$rsb$_',
    '{': '_$lcb$_',
    '}': '_$rcb$_',
}

HAIL_TYPE_TO_ES_TYPE_MAPPING = {
    hl.tint: "integer",
    hl.tint32: "integer",
    hl.tint64: "long",
    hl.tfloat: "double",
    hl.tfloat32: "float",
    hl.tfloat64: "double",
    hl.tstr: "keyword",
    hl.tbool: "boolean",
}

# used in wait_for_shard_transfer
# stolen from https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_client_v7.py#L20
LOADING_NODES_NAME = 'elasticsearch-es-data-loading*'


# https://hail.is/docs/devel/types.html
# https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-types.html
def _elasticsearch_mapping_for_type(dtype):
    """
    https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_utils.py#L53
    """
    if isinstance(dtype, hl.tstruct):
        return {"properties": {field: _elasticsearch_mapping_for_type(dtype[field]) for field in dtype.fields}}
    if isinstance(dtype, (hl.tarray, hl.tset)):
        element_mapping = _elasticsearch_mapping_for_type(dtype.element_type)
        if isinstance(dtype.element_type, hl.tstruct):
            element_mapping["type"] = "nested"
        return element_mapping
    if isinstance(dtype, hl.tlocus):
        return {"type": "object", "properties": {"contig": {"type": "keyword"}, "position": {"type": "integer"}}}
    if dtype in HAIL_TYPE_TO_ES_TYPE_MAPPING:
        return {"type": HAIL_TYPE_TO_ES_TYPE_MAPPING[dtype]}

    # tdict, ttuple, tlocus, tinterval, tcall
    raise NotImplementedError


def encode_field_name(s):
    """
    Encodes arbitrary string into an elasticsearch field name
    https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_utils.py#L123

    See:
    https://discuss.elastic.co/t/special-characters-in-field-names/10658/2
    https://discuss.elastic.co/t/illegal-characters-in-elasticsearch-field-names/17196/2
    """
    field_name = StringIO()
    for i, c in enumerate(s):
        if c == ES_FIELD_NAME_ESCAPE_CHAR:
            field_name.write(2 * ES_FIELD_NAME_ESCAPE_CHAR)
        elif c in ES_FIELD_NAME_SPECIAL_CHAR_MAP:
            field_name.write(ES_FIELD_NAME_SPECIAL_CHAR_MAP[c])  # encode the char
        else:
            field_name.write(c)  # write out the char as is

    field_name = field_name.getvalue()

    # escape 1st char if necessary
    if any(field_name.startswith(c) for c in ES_FIELD_NAME_BAD_LEADING_CHARS):
        return ES_FIELD_NAME_ESCAPE_CHAR + field_name
    else:
        return field_name


def struct_to_dict(struct):
    """
    https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/hail_elasticsearch_client.py#L21
    """
    return {k: dict(struct_to_dict(v)) if isinstance(v, hl.utils.Struct) else v for k, v in struct.items()}


class ElasticsearchClient:
    """
    Clone of https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/hail_elasticsearch_client.py#L25
    """

    def __init__(self, host: str, port: str, es_username: str, es_password: str):
        """
        This is a complete stripping of the ES Client in Hail Batch, which is thin-ish wrapper around a pile of
        hail methods
        """
        self._host = host
        self._port = port
        self._es_username = es_username
        self._es_password = es_password

        if not host.startswith('http://') or not host.startswith('https://'):
            host = f'https://{host}'
        _host = f'{host}:{port}'
        self.es = elasticsearch.Elasticsearch(_host, basic_auth=(self._es_username, self._es_password))

        # check connection
        logging.info(self.es.info())

    def wait_for_shard_transfer(self, index_name, num_attempts=1000):
        """
        Wait for shards to move off of the loading nodes before connecting to seqr
        https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_client_v7.py#L134
        """
        for i in range(num_attempts):
            shards = self.es.cat.shards(index=index_name)
            if LOADING_NODES_NAME not in shards:
                logging.warning("Shards are on {}".format(shards))
                return
            logging.warning(
                "Waiting for {} shards to transfer off the es-data-loading nodes: \n{}".format(
                    len(shards.strip().split("\n")),
                    shards,
                ),
            )
            time.sleep(5)

        raise Exception('Shards did not transfer off loading nodes')

    def create_mapping(self, index_name, elasticsearch_schema, num_shards=1, _meta=None):
        """
        Calls es.indices.create to create an elasticsearch index mapping.
        I've removed the alternate code paths for 'update' - we don't use it
        https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/hail_scripts/elasticsearch/elasticsearch_client_v7.py#L62

        Args:
            index_name (str): elasticsearch index mapping
            elasticsearch_schema (dict): elasticsearch mapping "properties" dictionary
            num_shards (int): how many shards the index will contain
            _meta (dict): optional _meta info for this index
                (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-meta-field.html)
        """

        index_mapping = {'properties': elasticsearch_schema}

        if _meta:
            logging.info(f'==> index _meta: {_meta}')
            index_mapping['_meta'] = _meta

        if not self.es.indices.exists(index=index_name):
            body = {
                'mappings': index_mapping,
                'settings': {
                    'number_of_shards': num_shards,
                    'number_of_replicas': 0,
                    'index.mapping.total_fields.limit': 10000,
                    'index.refresh_interval': -1,
                    'index.codec': 'best_compression',  # halves disk usage, no difference in query times
                },
            }

            logging.info(f'create_mapping - elasticsearch schema: \n{elasticsearch_schema}')
            logging.info('==> creating elasticsearch index {}'.format(index_name))

            self.es.indices.create(index=index_name, body=body)

    def export_table_to_elasticsearch(self, table, **kwargs):
        es_config = kwargs.get('elasticsearch_config', {})
        # to remove the write-null-values behaviour, remove the es.spark.dataframe.write.null entry
        es_config.update(
            {
                'es.net.ssl': 'true',
                'es.nodes.wan.only': 'true',
                'es.net.http.auth.user': self._es_username,
                'es.net.http.auth.pass': self._es_password,
                # investigate whether we lose anything from always writing Nulls
                # this is used in seqr-loading-pipelines, but not writing these would
                # result in a much smaller index
                'es.spark.dataframe.write.null': 'true',
                # We are not explicitly indexing the ES Index on varianId at this time
                # we should probably investigate this in future, but if we index on variantId
                # we run into the possibility that gCNV (currently multiple separate indices)
                # could have an ID clash, so the variant rows could overwrite each other
                # 'es.mapping.id': 'variantId'  # uncomment to explicitly index rows on the UID
            },
        )
        es_config['es.write.operation'] = 'index'
        # encode any special chars in column names
        rename_dict = {}
        for field_name in table.row_value.dtype.fields:
            encoded_name = field_name

            # optionally replace . with _ in a non-reversible way
            encoded_name = encoded_name.replace(".", '_')

            # replace all other special chars with an encoding that's uglier, but reversible
            encoded_name = encode_field_name(encoded_name)

            if encoded_name != field_name:
                rename_dict[field_name] = encoded_name

        for original_name, encoded_name in rename_dict.items():
            logging.info(f'Encoding column name {original_name} to {encoded_name}')

        table = table.rename(rename_dict)

        # create elasticsearch index with fields that match the ones in the table
        elasticsearch_schema = _elasticsearch_mapping_for_type(table.key_by().row_value.dtype)["properties"]

        index_name = kwargs['index_name']
        assert index_name

        if self.es.indices.exists(index=index_name):
            self.es.indices.delete(index=index_name)

        _meta = struct_to_dict(hl.eval(table.globals))

        self.create_mapping(index_name, elasticsearch_schema, num_shards=kwargs['num_shards'], _meta=_meta)

        hl.export_elasticsearch(table, self._host, int(self._port), index_name, '', 5000, es_config)
        self.es.indices.forcemerge(index=index_name, request_timeout=60)


def main():

    parser = ArgumentParser(description='Argument Parser for the ES generation script')
    parser.add_argument('--mt_path', help='MT path name', required=True)
    parser.add_argument('--index', help='ES index name', required=True)
    parser.add_argument('--flag', help='ES index "DONE" file path')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    password: str | None = read_secret(
        project_id=config_retrieve(['elasticsearch', 'password_project_id'], ''),
        secret_name=config_retrieve(['elasticsearch', 'password_secret_id'], ''),
        fail_gracefully=True,
    )

    # no password, but we fail gracefully
    if password is None:
        logging.warning(f'No permission to access ES password, skipping creation of {args.index}')
        exit(0)

    host = config_retrieve(['elasticsearch', 'host'])
    port = config_retrieve(['elasticsearch', 'port'])
    username = config_retrieve(['elasticsearch', 'username'])
    logging.info(f'Connecting to ElasticSearch: host="{host}", port="{port}", user="{username}"')

    ncpu = config_retrieve(['workflow', 'ncpu'], 4)
    hl.context.init_spark(master=f'local[{ncpu}]', quiet=True)
    hl.default_reference('GRCh38')

    mt = hl.read_matrix_table(args.mt_path)

    logging.info('Getting rows and exporting to the ES')

    # get the rows, flattened, stripped of key and VEP annotations
    row_ht = elasticsearch_row(mt)

    # Calculate the number of shards from the number of variants and sequencing groups.
    # https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/luigi_pipeline/lib/hail_tasks.py#L273
    # the denominator in this calculation used to be  1.4 * 10 ** 9, resulting in ~65GB shards
    # it's been reduced to give us more shards, closer to the optimum range 10-50GB
    # there's a huge inflation (3-4x) between the MT size on disk and the ES Index size, maybe this wasn't factored in?
    es_shards = math.ceil((mt.count_rows() * mt.count_cols()) / (10**9 / 2))

    es_client = ElasticsearchClient(host=host, port=port, es_username=username, es_password=password)

    # delete the index if it exists already - we shouldn't be doing this
    if es_client.es.indices.exists(index=args.index):
        es_client.es.indices.delete(index=args.index)

    es_client.export_table_to_elasticsearch(row_ht, index_name=args.index, num_shards=es_shards)

    # https://github.com/broadinstitute/seqr-loading-pipelines/blob/c113106204165e22b7a8c629054e94533615e7d2/luigi_pipeline/lib/hail_tasks.py#L266
    if es_shards < 25:
        es_client.wait_for_shard_transfer(args.index)

    with to_path(args.flag).open('w') as f:
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
    table = mt.rows()
    if 'vep' in table.row:
        table = table.drop('vep')
    key = table.key
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    flat_table = table.flatten()
    # When flattening, the table is unkeyed, which causes problems because our row keys should not
    # be normal fields. We can also re-key, but I believe this is computational?
    # PS: row key is often locus and allele, but does not have to be
    flat_table = flat_table.drop(*key)
    flat_table.describe()
    return flat_table


if __name__ == '__main__':
    main()
