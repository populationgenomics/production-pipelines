"""
Test Hail Query functions.
"""
import os
import string
import time
from random import choices
import logging
from os.path import exists

import toml
from cpg_utils import to_path, Path
from pytest_mock import MockFixture

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


def update_dict(d1: dict, d2: dict) -> None:
    """Updates the d1 dict with the values from the d2 dict recursively in-place."""
    for k, v2 in d2.items():
        v1 = d1.get(k)
        if isinstance(v1, dict) and isinstance(v2, dict):
            update_dict(v1, v2)
        else:
            d1[k] = v2


def timestamp(rand_suffix_len: int = 5) -> str:
    """
    Generate a timestamp string. If `rand_suffix_len` is set, adds a short random
    string of this length for uniqueness.
    """
    result = time.strftime('%Y_%m%d_%H%M')
    if rand_suffix_len:
        rand_bit = ''.join(
            choices(string.ascii_uppercase + string.digits, k=rand_suffix_len)
        )
        result += f'_{rand_bit}'
    return result


def _make_config(results_prefix: Path) -> dict:
    d: dict = {}
    for fp in (
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'workflows.toml',
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml',
    ):
        with fp.open():
            update_dict(d, toml.load(fp))

    update_dict(
        d,
        {
            'workflow': {
                'local_dir': str(results_prefix),
                'dataset_gcp_project': 'thousand-genomes',
                'dataset': 'thousand-genomes',
                'access_level': 'test',
                'sequencing_type': 'genome',
                'check_intermediates': True,
                'path_scheme': 'local',
                'reference_prefix': str(
                    to_path(__file__).parent / 'data' / 'large_cohort' / 'reference'
                ),
            },
            'large_cohort': {
                'sample_qc_cutoffs': {
                    'min_n_snps': 2500,  # to make it pass for toy subset
                },
                'n_pcs': 3,  # minimal allowed number
            },
            'hail': {
                'billing_project': 'thousand-genomes',
                'query_backend': 'spark_local',
            },
            'combiner': {
                'intervals': ['chr20:start-end', 'chrX:start-end', 'chrY:start-end'],
            },
        },
    )
    return d


def test_large_cohort(mocker: MockFixture):
    """
    Run entire workflow in a local mode.
    """
    results_prefix = (
        to_path(__file__).parent / 'results' / os.getenv('TEST_TIMESTAMP', timestamp())
    ).absolute()
    results_prefix.mkdir(parents=True, exist_ok=True)
    conf = _make_config(results_prefix)

    mocker.patch('cpg_utils.config.get_config', lambda: conf)

    from cpg_workflows.filetypes import GvcfPath
    from cpg_workflows.targets import Cohort

    cohort = Cohort()
    ds = cohort.create_dataset('thousand-genomes')
    gvcf_root = to_path(__file__).parent / 'data' / 'large_cohort' / 'gvcf'
    found_gvcf_paths = list(gvcf_root.glob('*.g.vcf.gz'))
    assert len(found_gvcf_paths) > 0, gvcf_root
    for gvcf_path in found_gvcf_paths:
        sample_id = gvcf_path.name.split('.')[0]
        s = ds.add_sample(id=sample_id, external_id=sample_id.replace('CPG', 'EXT'))
        s.gvcf = GvcfPath(gvcf_path)

    mocker.patch('cpg_workflows.inputs.create_cohort', lambda: cohort)

    from cpg_workflows.large_cohort import combiner
    from cpg_workflows.large_cohort import ancestry_pca
    from cpg_workflows.large_cohort import sample_qc
    from cpg_workflows.large_cohort import dense_subset
    from cpg_workflows.large_cohort import relatedness
    from cpg_workflows.large_cohort import ancestry_plots
    from cpg_workflows.large_cohort import site_only_vcf
    from cpg_utils.hail_batch import start_query_context

    start_query_context()

    vds_path = results_prefix / 'v01.vds'
    combiner.run(out_vds_path=vds_path, tmp_prefix=results_prefix / 'tmp')

    sample_qc_ht_path = results_prefix / 'sample_qc.ht'
    sample_qc.run(
        vds_path=vds_path,
        out_sample_qc_ht_path=sample_qc_ht_path,
        tmp_prefix=results_prefix / 'tmp',
    )

    dense_mt_path = results_prefix / 'dense.mt'
    dense_subset.run(
        vds_path=vds_path,
        out_dense_mt_path=dense_mt_path,
    )

    relateds_to_drop_ht_path = results_prefix / 'relateds_to_drop.ht'
    relatedness.run(
        dense_mt_path=dense_mt_path,
        sample_qc_ht_path=sample_qc_ht_path,
        out_relatedness_ht_path=results_prefix / 'relatedness.ht',
        out_relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=results_prefix / 'tmp',
    )

    scores_ht_path = results_prefix / 'scores.ht'
    eigenvalues_ht_path = results_prefix / 'eigenvalues.ht'
    loadings_ht_path = results_prefix / 'loadings.ht'
    inferred_pop_ht_path = results_prefix / 'inferred_pop.ht'
    ancestry_pca.run(
        dense_mt_path=dense_mt_path,
        sample_qc_ht_path=sample_qc_ht_path,
        relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=results_prefix / 'tmp',
        out_scores_ht_path=scores_ht_path,
        out_eigenvalues_ht_path=eigenvalues_ht_path,
        out_loadings_ht_path=loadings_ht_path,
        out_inferred_pop_ht_path=inferred_pop_ht_path,
    )
    ancestry_plots.run(
        out_path_pattern=results_prefix / 'plots' / '{scope}_pc{pci}.{ext}',
        sample_qc_ht_path=sample_qc_ht_path,
        scores_ht_path=scores_ht_path,
        eigenvalues_ht_path=eigenvalues_ht_path,
        loadings_ht_path=loadings_ht_path,
        inferred_pop_ht_path=inferred_pop_ht_path,
    )

    siteonly_vcf_path = results_prefix / 'siteonly.vcf.bgz'
    site_only_vcf.run(
        vds_path=vds_path,
        sample_qc_ht_path=sample_qc_ht_path,
        relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        out_vcf_path=siteonly_vcf_path,
        tmp_prefix=results_prefix / 'tmp',
    )

    assert exists(vds_path)
    assert exists(results_prefix / 'plots' / 'dataset_pc1.html')
    assert exists(siteonly_vcf_path)
