"""
Utility functions for Query script submitted with Dataproc. Do not depend on
cpg-pipes and python 3.10.
"""

import asyncio
import os
import typing

import hail as hl
from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build, dataset_path


def start_hail_context(
    query_backend: typing.Literal['spark', 'batch', 'local', 'spark_local']
    | None = None,
    log_path: str | None = None,
    dataset: str | None = None,
    billing_project: str | None = None,
):
    """
    Start Hail context, depending on the Backend class in hail/query_backend
    """
    query_backend = query_backend or get_config().get('hail', {}).get(
        'query_backend', 'spark'
    )
    if query_backend == 'spark':
        hl.init(default_reference=genome_build())
    elif query_backend == 'spark_local':
        hl.init(
            default_reference=genome_build(),
            master='local[2]',
            quiet=True,
            log=log_path or dataset_path('hail-log.txt', category='tmp'),
        )
    elif query_backend == 'local':
        hl.utils.java.Env.hc()  # force initialization
    else:
        assert query_backend == 'batch'
        if hl.utils.java.Env._hc:  # pylint: disable=W0212
            return  # already initialised
        dataset = dataset or get_config()['workflow']['dataset']
        billing_project = billing_project or get_config()['hail']['billing_project']

        asyncio.get_event_loop().run_until_complete(
            hl.init_batch(
                billing_project=billing_project,
                remote_tmpdir=f'gs://cpg-{dataset}-hail/batch-tmp',
                token=os.environ.get('HAIL_TOKEN'),
                default_reference='GRCh38',
            )
        )


# def get_vds(
#     vds_path: str,
#     split: bool = False,
#     sample_ht: Optional[hl.Table] = None,
#     pass_only: bool = False,
#     n_partitions: int = None,
# ) -> hl.vds.VariantDataset:
#     """
#     Wrapper function to get VDS with desired filtering and metadata annotations.
#     @param vds_path: path to VDS
#     @param split:
#         Split multiallelics and convert local-allele LGT/LA fields to GT.
#         Note: this will perform a split on the MT rather than grab an already split MT
#     @param sample_ht: sample-level metadata. Required for pass_only and release_only
#         filters.
#     @param pass_only: remove samples that failed filtering requires sample_ht.
#     @param n_partitions: number of partitions to use to load the VDS.
#     @return: VDS with chosen annotations and filters.
#     """
#     vds = hl.vds.read_vds(vds_path, n_partitions=n_partitions)
#
#     if pass_only:
#         assert sample_ht is not None
#         failed_sample_ht = sample_ht.filter(hl.len(sample_ht.filters) > 0)
#         vds = hl.vds.filter_samples(
#             vds,
#             failed_sample_ht,
#             keep=False,
#             remove_dead_alleles=True,
#         )
#
#     if split:
#         vmt = vds.variant_data
#         vmt = vmt.annotate_rows(
#             n_unsplit_alleles=hl.len(vmt.alleles),
#             mixed_site=(hl.len(vmt.alleles) > 2)
#             & hl.any(lambda a: hl.is_indel(vmt.alleles[0], a), vmt.alleles[1:])
#             & hl.any(lambda a: hl.is_snp(vmt.alleles[0], a), vmt.alleles[1:]),
#         )
#         vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
#         vds = hl.vds.VariantDataset(vds.reference_data, vmt)
#
#     return vds
