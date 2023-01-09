"""
Test large-cohort workflow.
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
from . import results_prefix, update_dict

ref_prefix = to_path(__file__).parent / 'data/large_cohort/reference'
gnomad_prefix = ref_prefix / 'gnomad/v0'
broad_prefix = ref_prefix / 'hg38/v0'

TOML = f"""
[workflow]
dataset_gcp_project = "thousand-genomes"
dataset = "thousand-genomes"
access_level = "test"
sequencing_type = "genome"
check_intermediates = true

[large_cohort]
n_pcs = 3

[hail]
billing_project = "thousand-genomes"
query_backend = "spark_local"

[combiner]
intervals = [ "chr20:start-end", "chrX:start-end", "chrY:start-end",]

[references]
genome_build = "GRCh38"

[storage.default]
default = "{results_prefix()}"
web = "{results_prefix()}-web"
analysis = "{results_prefix()}-analysis"
tmp = "{results_prefix()}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.thousand-genomes]
default = "{results_prefix()}"
web = "{results_prefix()}-web"
analysis = "{results_prefix()}-analysis"
tmp = "{results_prefix()}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[large_cohort.sample_qc_cutoffs]
min_n_snps = 2500

[references.gnomad]
tel_and_cent_ht = "{gnomad_prefix}/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht"
lcr_intervals_ht = "{gnomad_prefix}/lcr_intervals/LCRFromHengHg38.ht"
seg_dup_intervals_ht = "{gnomad_prefix}/seg_dup_intervals/GRCh38_segdups.ht"
clinvar_ht = "{gnomad_prefix}/clinvar/clinvar_20190923.ht"
hapmap_ht = "{gnomad_prefix}/hapmap/hapmap_3.3.hg38.ht"
kgp_omni_ht = "{gnomad_prefix}/kgp/1000G_omni2.5.hg38.ht"
kgp_hc_ht = "{gnomad_prefix}/kgp/1000G_phase1.snps.high_confidence.hg38.ht"
mills_ht = "{gnomad_prefix}/mills/Mills_and_1000G_gold_standard.indels.hg38.ht"
predetermined_qc_variants = "{gnomad_prefix}/sample_qc/pre_ld_pruning_qc_variants.ht"

[references.broad]
genome_calling_interval_lists = "{broad_prefix}/wgs_calling_regions.hg38.interval_list"
protein_coding_gtf = "{broad_prefix}/sv-resources/resources/v1/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf"
"""


def _mock_config() -> dict:
    d: dict = {}
    for fp in [
        to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml',
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml',
    ]:
        with fp.open():
            update_dict(d, toml.load(fp))

    update_dict(d, toml.loads(TOML))
    return d


def _mock_cohort():
    from cpg_workflows.filetypes import GvcfPath
    from cpg_workflows.targets import Cohort

    cohort = Cohort()
    dataset = cohort.create_dataset('thousand-genomes')
    gvcf_root = to_path(__file__).parent / 'data' / 'large_cohort' / 'gvcf'
    found_gvcf_paths = list(gvcf_root.glob('*.g.vcf.gz'))
    assert len(found_gvcf_paths) > 0, gvcf_root
    for gvcf_path in found_gvcf_paths:
        sample_id = gvcf_path.name.split('.')[0]
        sample = dataset.add_sample(
            id=sample_id, external_id=sample_id.replace('CPG', 'EXT')
        )
        sample.gvcf = GvcfPath(gvcf_path)
    return cohort


def test_large_cohort(mocker: MockFixture):
    """
    Run entire workflow in a local mode.
    """
    mocker.patch('cpg_utils.config.get_config', _mock_config)
    mocker.patch('cpg_workflows.inputs.create_cohort', _mock_cohort)

    from cpg_workflows.large_cohort import combiner
    from cpg_workflows.large_cohort import ancestry_pca
    from cpg_workflows.large_cohort import sample_qc
    from cpg_workflows.large_cohort import dense_subset
    from cpg_workflows.large_cohort import relatedness
    from cpg_workflows.large_cohort import ancestry_plots
    from cpg_workflows.large_cohort import site_only_vcf
    from cpg_utils.hail_batch import start_query_context

    start_query_context()

    res_pref = to_path(results_prefix())
    vds_path = res_pref / 'v01.vds'
    combiner.run(out_vds_path=vds_path, tmp_prefix=res_pref / 'tmp')

    sample_qc_ht_path = res_pref / 'sample_qc.ht'
    sample_qc.run(
        vds_path=vds_path,
        out_sample_qc_ht_path=sample_qc_ht_path,
        tmp_prefix=res_pref / 'tmp',
    )

    dense_mt_path = res_pref / 'dense.mt'
    dense_subset.run(
        vds_path=vds_path,
        out_dense_mt_path=dense_mt_path,
    )

    relateds_to_drop_ht_path = res_pref / 'relateds_to_drop.ht'
    relatedness.run(
        dense_mt_path=dense_mt_path,
        sample_qc_ht_path=sample_qc_ht_path,
        out_relatedness_ht_path=res_pref / 'relatedness.ht',
        out_relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=res_pref / 'tmp',
    )

    scores_ht_path = res_pref / 'scores.ht'
    eigenvalues_ht_path = res_pref / 'eigenvalues.ht'
    loadings_ht_path = res_pref / 'loadings.ht'
    inferred_pop_ht_path = res_pref / 'inferred_pop.ht'
    ancestry_pca.run(
        dense_mt_path=dense_mt_path,
        sample_qc_ht_path=sample_qc_ht_path,
        relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=res_pref / 'tmp',
        out_scores_ht_path=scores_ht_path,
        out_eigenvalues_ht_path=eigenvalues_ht_path,
        out_loadings_ht_path=loadings_ht_path,
        out_inferred_pop_ht_path=inferred_pop_ht_path,
    )
    ancestry_plots.run(
        out_path_pattern=res_pref / 'plots' / '{scope}_pc{pci}.{ext}',
        sample_qc_ht_path=sample_qc_ht_path,
        scores_ht_path=scores_ht_path,
        eigenvalues_ht_path=eigenvalues_ht_path,
        loadings_ht_path=loadings_ht_path,
        inferred_pop_ht_path=inferred_pop_ht_path,
    )

    siteonly_vcf_path = res_pref / 'siteonly.vcf.bgz'
    site_only_vcf.run(
        vds_path=vds_path,
        sample_qc_ht_path=sample_qc_ht_path,
        relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        out_vcf_path=siteonly_vcf_path,
        tmp_prefix=res_pref / 'tmp',
    )

    assert exists(vds_path)
    assert exists(res_pref / 'plots' / 'dataset_pc1.html')
    assert exists(siteonly_vcf_path)
