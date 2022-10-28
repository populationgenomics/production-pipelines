"""
Test Hail Query functions.
"""
import os
import logging
from os.path import exists

import toml
from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths
from cpg_utils.config import update_dict
from pytest_mock import MockFixture

from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import timestamp

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


def _set_config(results_prefix: Path, extra_conf: dict | None = None):
    d = {
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
            'dry_run': True,
            'query_backend': 'spark_local',
        },
    }
    if extra_conf:
        update_dict(d, extra_conf)
    with (conf_path := results_prefix / 'config.toml').open('w') as f:
        toml.dump(d, f)

    set_config_paths(
        [
            str(p)
            for p in [
                to_path(__file__).parent.parent
                / 'configs'
                / 'defaults'
                / 'workflows.toml',
                to_path(__file__).parent.parent
                / 'configs'
                / 'defaults'
                / 'large_cohort.toml',
                conf_path,
            ]
        ]
    )


def test_large_cohort(mocker: MockFixture):
    """
    Run entire workflow in a local mode.
    """
    results_prefix = (
        to_path(__file__).parent / 'results' / os.getenv('TEST_TIMESTAMP', timestamp())
    ).absolute()
    results_prefix.mkdir(parents=True, exist_ok=True)

    _set_config(
        results_prefix=results_prefix,
        extra_conf={
            'combiner': {
                'intervals': ['chr20:start-end', 'chrX:start-end', 'chrY:start-end'],
            },
        },
    )

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
