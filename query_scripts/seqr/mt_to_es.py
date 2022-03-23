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
        index_type_name='',  # need to be empty for ES=8.x.x
    )

    mt = hl.read_matrix_table(mt_path)
    
    # AS_MQ and InbreedingCoeff can be NaN, need to fix that to avoid failures like 
    # `is.hail.relocated.org.elasticsearch.hadoop.rest.EsHadoopRemoteException: 
    # mapper_parsing_exception: failed to parse field [info_AS_MQ] of type [double]`
    # in: https://batch.hail.populationgenomics.org.au/batches/6621/jobs/6
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AS_MQ=hl.if_else(
                hl.is_nan(mt.info.AS_MQ), 
                0.0, 
                mt.info.AS_MQ
            ),
            InbreedingCoeff=hl.if_else(
                hl.is_nan(mt.info.InbreedingCoeff), 
                0.0, 
                mt.info.InbreedingCoeff
            ),
            # https://batch.hail.populationgenomics.org.au/batches/6973/jobs/12
            AS_InbreedingCoeff=hl.if_else(
                hl.is_nan(mt.info.AS_InbreedingCoeff), 
                0.0, 
                mt.info.AS_InbreedingCoeff
            ),
        )
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
    """
    # Converts a mt to the row equivalent.
    ht = mt.rows()
    # Converts nested structs into one field, e.g. {a: {b: 1}} => a.b: 1
    table = ht.drop('vep').flatten()
    # When flattening, the table is unkeyed, which causes problems because our locus and alleles should not
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






# import re
# import logging
# import time
# from pprint import pformat
# 
# try:
#     import elasticsearch
# except ImportError:
#     import os
#     os.system("pip install elasticsearch>=8,<9")
#     import elasticsearch
# 
# handlers = set(logging.root.handlers)
# logging.root.handlers = list(handlers)
# logger = logging.getLogger()
# 
# LOADING_NODES_NAME = 'elasticsearch-es-data-loading*'
# 
# 
# class ElasticsearchClient:
# 
#     def __init__(self, host='localhost', port='9200', es_username='pipeline',
#                  es_password=None, es_use_ssl=False):
#         """Constructor.
# 
#         Args:
#             host (str): Elasticsearch server host
#             port (str): Elasticsearch server port
#             es_username (str): Elasticsearch username
#             es_password (str): Elasticsearch password
#             es_use_ssl (bool): Whether to use SSL for ElasticSearch connections
#         """
# 
#         self._es_username = es_username
#         self._es_password = es_password
#         self._es_use_ssl = es_use_ssl
#         self._host = host
#         self._port = port
# 
#         auth = (self._es_username, self._es_password) if self._es_password else None
# 
#         if not host.startswith('http://') or not host.startswith('https://'):
#             scheme = 'https' if es_use_ssl else 'http'
#             host = f'{scheme}://{host}'
#         _host = f'{host}:{port}'
#         self.es = elasticsearch.Elasticsearch(_host, basic_auth=auth)
# 
#         # check connection
#         logger.info(pformat(self.es.info()))
# 
#     def create_index(self, index_name, elasticsearch_schema, num_shards=1, _meta=None):
#         """Calls es.indices.create to create an elasticsearch index with the 
#         appropriate mapping.
# 
#         Args:
#             index_name (str): elasticsearch index mapping
#             elasticsearch_schema (dict): elasticsearch mapping "properties" dictionary
#             num_shards (int): how many shards the index will contain
#             _meta (dict): optional _meta info for this index
#                 (see https://www.elastic.co/guide/en/elasticsearch/reference/current
#                 /mapping-meta-field.html)
#         """
# 
#         self.create_or_update_mapping(
#             index_name, elasticsearch_schema, num_shards=num_shards, _meta=_meta,
#             create_only=True,
#         )
# 
#     def create_or_update_mapping(self, index_name, elasticsearch_schema, num_shards=1,
#                                  _meta=None, create_only=False):
#         """Calls es.indices.create or es.indices.put_mapping to create or update an 
#         elasticsearch index mapping.
# 
#         Args:
#             index_name (str): elasticsearch index mapping
#             elasticsearch_schema (dict): elasticsearch mapping "properties" dictionary
#             num_shards (int): how many shards the index will contain
#             _meta (dict): optional _meta info for this index
#                 (see https://www.elastic.co/guide/en/elasticsearch/reference/current
#                 /mapping-meta-field.html)
#             create_only (bool): only allow index creation, throws an error if index 
#             already exists
#         """
# 
#         index_mapping = {
#             'properties': elasticsearch_schema,
#         }
# 
#         if _meta:
#             logger.info('==> index _meta: ' + pformat(_meta))
#             index_mapping['_meta'] = _meta
# 
#         if not self.es.indices.exists(index=index_name):
#             body = {
#                 'mappings': index_mapping,
#                 'settings': {
#                     'number_of_shards': num_shards,
#                     'number_of_replicas': 0,
#                     'index.mapping.total_fields.limit': 10000,
#                     'index.refresh_interval': -1,
#                     'index.codec': 'best_compression',
#                     # halves disk usage, no difference in query times
#                 }
#             }
# 
#             logger.info('create_mapping - elasticsearch schema: \n' + pformat(
#                 elasticsearch_schema))
#             logger.info('==> creating elasticsearch index {}'.format(index_name))
# 
#             self.es.indices.create(index=index_name, body=body)
#         else:
#             if create_only:
#                 raise ValueError('Index {} already exists'.format(index_name))
# 
#             logger.info('==> updating elasticsearch index {}. New schema:\n{}'.format(
#                 index_name, pformat(elasticsearch_schema)))
# 
#             self.es.indices.put_mapping(index=index_name, body=index_mapping)
# 
#     def route_index_to_temp_es_cluster(self, index_name):
#         """
#         Apply shard allocation filtering rules to route the given index to 
#         elasticsearch loading nodes
#         """
#         self._update_settings(index_name, {
#             'index.routing.allocation.require._name': LOADING_NODES_NAME,
#             'index.routing.allocation.exclude._name': ''
#         })
# 
#     def route_index_off_temp_es_cluster(self, index_name):
#         """
#         Move any shards in the given index off of loading nodes
#         """
#         self._update_settings(index_name, {
#             'index.routing.allocation.require._name': '',
#             'index.routing.allocation.exclude._name': LOADING_NODES_NAME
#         })
# 
#     def _update_settings(self, index_name, body):
#         logger.info('==> Setting {} settings = {}'.format(index_name, body))
# 
#         self.es.indices.put_settings(index=index_name, body=body)
# 
#     def get_index_meta(self, index_name):
#         mappings = self.es.indices.get_mapping(index=index_name)
#         return mappings.get(index_name, {}).get('mappings', {}).get('_meta', {})
# 
#     def wait_for_shard_transfer(self, index_name, num_attempts=1000):
#         """
#         Wait for shards to move off of the loading nodes before connecting to seqr 
#         """
#         for i in range(num_attempts):
#             shards = self.es.cat.shards(index=index_name)
#             if LOADING_NODES_NAME not in shards:
#                 logger.warning("Shards are on {}".format(shards))
#                 return
#             logger.warning(
#                 "Waiting for {} shards to transfer off the es-data-loading nodes: \n{"
#                 "}".format(
#                     len(shards.strip().split("\n")), shards))
#             time.sleep(5)
# 
#         raise Exception('Shards did not transfer off loading nodes')
# 
# from io import StringIO
# 
# logger = logging.getLogger()
# 
# # Elastic search write operations.
# # See https://www.elastic.co/guide/en/elasticsearch/hadoop/current/configuration.html#_operation
# ELASTICSEARCH_INDEX = "index"
# ELASTICSEARCH_CREATE = "create"
# ELASTICSEARCH_UPDATE = "update"
# ELASTICSEARCH_UPSERT = "upsert"
# ELASTICSEARCH_WRITE_OPERATIONS = set([
#     ELASTICSEARCH_INDEX, ELASTICSEARCH_CREATE, ELASTICSEARCH_UPDATE, ELASTICSEARCH_UPSERT,
# ])
# 
# # make encoded values as human-readable as possible
# ES_FIELD_NAME_ESCAPE_CHAR = '$'
# ES_FIELD_NAME_BAD_LEADING_CHARS = set(['_', '-', '+', ES_FIELD_NAME_ESCAPE_CHAR])
# ES_FIELD_NAME_SPECIAL_CHAR_MAP = {
#     '.': '_$dot$_',
#     ',': '_$comma$_',
#     '#': '_$hash$_',
#     '*': '_$star$_',
#     '(': '_$lp$_',
#     ')': '_$rp$_',
#     '[': '_$lsb$_',
#     ']': '_$rsb$_',
#     '{': '_$lcb$_',
#     '}': '_$rcb$_',
# }
# 
# HAIL_TYPE_TO_ES_TYPE_MAPPING = {
#     hl.tint: "integer",
#     hl.tint32: "integer",
#     hl.tint64: "long",
#     hl.tfloat: "double",
#     hl.tfloat32: "float",
#     hl.tfloat64: "double",
#     hl.tstr: "keyword",
#     hl.tbool: "boolean",
# }
# 
# 
# # https://hail.is/docs/devel/types.html
# # https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-types.html
# def _elasticsearch_mapping_for_type(dtype):
#     if isinstance(dtype, hl.tstruct):
#         return {"properties": {field: _elasticsearch_mapping_for_type(dtype[field]) for field in dtype.fields}}
#     if isinstance(dtype, (hl.tarray, hl.tset)):
#         element_mapping = _elasticsearch_mapping_for_type(dtype.element_type)
#         if isinstance(dtype.element_type, hl.tstruct):
#             element_mapping["type"] = "nested"
#         return element_mapping
#     if isinstance(dtype, hl.tlocus):
#         return {"type": "object", "properties": {"contig": {"type": "keyword"}, "position": {"type": "integer"}}}
#     if dtype in HAIL_TYPE_TO_ES_TYPE_MAPPING:
#         return {"type": HAIL_TYPE_TO_ES_TYPE_MAPPING[dtype]}
# 
#     # tdict, ttuple, tlocus, tinterval, tcall
#     raise NotImplementedError
# 
# 
# def elasticsearch_schema_for_table(table, disable_doc_values_for_fields=(), disable_index_for_fields=()):
#     """
#     Converts the type of a table's row values into a dictionary that can be plugged in to
#     an elasticsearch mapping definition.
# 
#     Args:
#         table (hail.Table): the table to generate a schema for
#         disable_doc_values_for_fields: (optional) list of field names (the way they will be
#             named in the elasticsearch index) for which to not store doc_values
#             (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-params.html)
#         disable_index_for_fields: (optional) list of field names (the way they will be
#             named in the elasticsearch index) that shouldn't be indexed
#             (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-params.html)
#     Returns:
#         A dict that can be plugged in to an elasticsearch mapping as the value for "properties".
#         (see https://www.elastic.co/guide/en/elasticsearch/guide/current/root-object.html)
#     """
#     properties = _elasticsearch_mapping_for_type(table.key_by().row_value.dtype)["properties"]
# 
#     if disable_doc_values_for_fields:
#         logger.info("==> will disable doc values for %s", ", ".join(disable_doc_values_for_fields))
#         for es_field_name in disable_doc_values_for_fields:
#             if es_field_name not in properties:
#                 raise ValueError(
#                     "'%s' in disable_doc_values_for_fields arg is not in the elasticsearch schema: %s"
#                     % (es_field_name, properties)
#                 )
#             properties[es_field_name]["doc_values"] = False
# 
#     if disable_index_for_fields:
#         logger.info("==> will disable index fields for %s", ", ".join(disable_index_for_fields))
#         for es_field_name in disable_index_for_fields:
#             if es_field_name not in properties:
#                 flattened_fields = [key for key in properties if key.startswith(f'{es_field_name}_')]
#                 if flattened_fields:
#                     for flattened_es_field_name in flattened_fields:
#                         properties[flattened_es_field_name]["index"] = False
#                 else:
#                     raise ValueError(
#                         "'%s' in disable_index_for_fields arg is not in the elasticsearch schema: %s"
#                         % (es_field_name, properties)
#                     )
#             else:
#                 properties[es_field_name]["index"] = False
# 
#     return properties
# 
# def encode_field_name(s):
#     """Encodes arbitrary string into an elasticsearch field name
# 
#     See:
#     https://discuss.elastic.co/t/special-characters-in-field-names/10658/2
#     https://discuss.elastic.co/t/illegal-characters-in-elasticsearch-field-names/17196/2
#     """
#     field_name = StringIO()
#     for i, c in enumerate(s):
#         if c == ES_FIELD_NAME_ESCAPE_CHAR:
#             field_name.write(2*ES_FIELD_NAME_ESCAPE_CHAR)
#         elif c in ES_FIELD_NAME_SPECIAL_CHAR_MAP:
#             field_name.write(ES_FIELD_NAME_SPECIAL_CHAR_MAP[c])  # encode the char
#         else:
#             field_name.write(c)  # write out the char as is
# 
#     field_name = field_name.getvalue()
# 
#     # escape 1st char if necessary
#     if any(field_name.startswith(c) for c in ES_FIELD_NAME_BAD_LEADING_CHARS):
#         return ES_FIELD_NAME_ESCAPE_CHAR + field_name
#     else:
#         return field_name
# 
# 
# 
# def struct_to_dict(struct):
#     return {k: dict(struct_to_dict(v)) if isinstance(v, hl.utils.Struct) else v for k, v in struct.items()}
# 
# 
# class HailElasticsearchClient(ElasticsearchClient):
#     def export_table_to_elasticsearch(
#         self,
#         table: hl.Table,
#         index_name :str = "data",
#         index_type_name :str = '_doc',
#         block_size :int = 5000,
#         num_shards :int = 10,
#         delete_index_before_exporting :bool = True,
#         elasticsearch_write_operation :str = ELASTICSEARCH_INDEX,
#         ignore_elasticsearch_write_errors :bool = False,
#         elasticsearch_mapping_id=None,
#         field_name_to_elasticsearch_type_map=None,
#         disable_doc_values_for_fields=(),
#         disable_index_for_fields=(),
#         field_names_replace_dot_with="_",
#         func_to_run_after_index_exists=None,
#         export_globals_to_index_meta=True,
#         verbose=True,
#         write_null_values=False,
#         elasticsearch_config=None,
#     ):
#         """Create a new elasticsearch index to store the records in this table, and then export all records to it.
# 
#         Args:
#             table (Table): hail Table
#             index_name (string): elasticsearch index name
#             index_type_name (string): elasticsearch index type
#             block_size (int): number of records to write in one bulk insert
#             num_shards (int): number of shards to use for this index
#                 (see https://www.elastic.co/guide/en/elasticsearch/guide/current/overallocation.html)
#             delete_index_before_exporting (bool): Whether to drop and re-create the index before exporting.
#             elasticsearch_write_operation (string): Can be one of these constants:
#                     ELASTICSEARCH_INDEX
#                     ELASTICSEARCH_CREATE
#                     ELASTICSEARCH_UPDATE
#                     ELASTICSEARCH_UPSERT
#                 See https://www.elastic.co/guide/en/elasticsearch/hadoop/current/configuration.html#_operation
#             ignore_elasticsearch_write_errors (bool): If True, elasticsearch errors will be logged, but won't cause
#                 the bulk write call to throw an error. This is useful when, for example,
#                 elasticsearch_write_operation="update", and the desired behavior is to update all documents that exist,
#                 but to ignore errors for documents that don't exist.
#             elasticsearch_mapping_id (str): if specified, sets the es.mapping.id which is the column name to use as the document ID
#                 See https://www.elastic.co/guide/en/elasticsearch/hadoop/current/configuration.html#cfg-mapping
#             field_name_to_elasticsearch_type_map (dict): (optional) a map of table field names to
#                 their elasticsearch field spec - for example: {
#                     'allele_freq': { 'type': 'half_float' },
#                     ...
#                 }.
#                 See https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping.html for
#                 more details. Any values in this dictionary will override
#                 the default type mapping derived from the hail table's row type.
#                 Field names can be regular expressions.
#             disable_doc_values_for_fields (tuple): (optional) list of field names (the way they will be
#                 named in the elasticsearch index) for which to not store doc_values
#                 (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-params.html)
#             disable_index_for_fields (tuple): (optional) list of field names (the way they will be
#                 named in the elasticsearch index) that shouldn't be indexed
#                 (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-params.html)
#             field_names_replace_dot_with (string): since "." chars in field names are interpreted in
#                 special ways by elasticsearch, set this arg to first go through and replace "." with
#                 this string in all field names. This replacement is not reversible (or atleast not
#                 unambiguously in the general case) Set this to None to disable replacement, and fall back
#                 on an encoding that's uglier, but reversible (eg. "." will be converted to "_$dot$_")
#             func_to_run_after_index_exists (function): optional function to run after creating the index, but before exporting any data.
#             export_globals_to_index_meta (bool): whether to add table.globals object to the index _meta field:
#                 (see https://www.elastic.co/guide/en/elasticsearch/reference/current/mapping-meta-field.html)
#             verbose (bool): whether to print schema and stats
#             write_null_values (bool): whether to write fields that are null to the index
#             elasticsearch_config: The initial elasticsearch config from the caller
#         """
# 
#         elasticsearch_config = elasticsearch_config or {}
#         if (
#             elasticsearch_write_operation is not None
#             and elasticsearch_write_operation not in ELASTICSEARCH_WRITE_OPERATIONS
#         ):
#             raise ValueError(
#                 "Unexpected value for elasticsearch_write_operation arg: " + str(elasticsearch_write_operation)
#             )
# 
#         if elasticsearch_write_operation is not None:
#             elasticsearch_config["es.write.operation"] = elasticsearch_write_operation
# 
#         if elasticsearch_write_operation in (ELASTICSEARCH_UPDATE, ELASTICSEARCH_UPSERT) or write_null_values:
#             # see https://www.elastic.co/guide/en/elasticsearch/hadoop/master/spark.html#spark-sql-write
#             # "By default, elasticsearch-hadoop will ignore null values in favor of not writing any field at all.
#             # If updating/upserting, then existing field values may need to be overwritten with nulls
#             elasticsearch_config["es.spark.dataframe.write.null"] = "true"
# 
#         if elasticsearch_mapping_id is not None:
#             elasticsearch_config["es.mapping.id"] = elasticsearch_mapping_id
# 
#         if ignore_elasticsearch_write_errors:
#             # see docs in https://www.elastic.co/guide/en/elasticsearch/hadoop/current/errorhandlers.html
#             elasticsearch_config["es.write.rest.error.handlers"] = "log"
#             elasticsearch_config["es.write.rest.error.handler.log.logger.name"] = "BulkErrors"
# 
#         if self._es_password:
#             elasticsearch_config.update({
#                 'es.net.http.auth.user': self._es_username,
#                 'es.net.http.auth.pass': self._es_password,
#             })
# 
#         if self._es_use_ssl:
#             elasticsearch_config['es.net.ssl'] = 'true'
#             # If using SSL, the instance is likely managed, in which case we
#             # can't discover nodes.
#             elasticsearch_config['es.nodes.wan.only'] = 'true'
# 
#         # encode any special chars in column names
#         rename_dict = {}
#         for field_name in table.row_value.dtype.fields:
#             encoded_name = field_name
# 
#             # optionally replace . with _ in a non-reversible way
#             if field_names_replace_dot_with is not None:
#                 encoded_name = encoded_name.replace(".", field_names_replace_dot_with)
# 
#             # replace all other special chars with an encoding that's uglier, but reversible
#             encoded_name = encode_field_name(encoded_name)
# 
#             if encoded_name != field_name:
#                 rename_dict[field_name] = encoded_name
# 
#         for original_name, encoded_name in rename_dict.items():
#             logger.info("Encoding column name %s to %s", original_name, encoded_name)
# 
#         table = table.rename(rename_dict)
# 
#         if verbose:
#             logger.info(pformat(table.row_value.dtype))
# 
#         # create elasticsearch index with fields that match the ones in the table
#         elasticsearch_schema = elasticsearch_schema_for_table(
#             table,
#             disable_doc_values_for_fields=disable_doc_values_for_fields,
#             disable_index_for_fields=disable_index_for_fields,
#         )
# 
#         # override elasticsearch types
#         if field_name_to_elasticsearch_type_map is not None:
#             modified_elasticsearch_schema = dict(elasticsearch_schema)  # make a copy
#             for field_name_regexp, elasticsearch_field_spec in field_name_to_elasticsearch_type_map.items():
#                 match_count = 0
#                 for key in elasticsearch_schema.keys():
#                     if re.match(field_name_regexp, key):
#                         modified_elasticsearch_schema[key] = elasticsearch_field_spec
#                         match_count += 1
# 
#                 logger.info("%d columns matched '%s'", match_count, field_name_regexp)
# 
#             elasticsearch_schema = modified_elasticsearch_schema
# 
#         # optionally delete the index before creating it
#         if delete_index_before_exporting and self.es.indices.exists(index=index_name):
#             self.es.indices.delete(index=index_name)
# 
#         _meta = None
#         if export_globals_to_index_meta:
#             _meta = struct_to_dict(hl.eval(table.globals))
# 
#         self.create_or_update_mapping(index_name, elasticsearch_schema, num_shards=num_shards, _meta=_meta)
# 
#         if func_to_run_after_index_exists:
#             func_to_run_after_index_exists()
# 
#         logger.info(
#             "==> exporting data to elasticsearch. Write mode: %s, blocksize: %d",
#             elasticsearch_write_operation,
#             block_size,
#         )
# 
#         hl.export_elasticsearch(
#             table, self._host, int(self._port), index_name, '', block_size, elasticsearch_config, verbose
#         )
# 
#         """
#         Potentially useful config settings for export_elasticsearch(..)
#         (https://www.elastic.co/guide/en/elasticsearch/hadoop/current/configuration.html)
# 
#         es.write.operation // default: index (create, update, upsert)
#         es.http.timeout // default 1m
#         es.http.retries // default 3
#         es.batch.size.bytes  // default 1mb
#         es.batch.size.entries  // default 1000
#         es.batch.write.refresh // default true  (Whether to invoke an index refresh or not after a bulk update has been completed)
#         """
# 
#         self.es.indices.forcemerge(index=index_name, request_timeout=60)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
