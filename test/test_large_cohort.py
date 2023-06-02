"""
Test large-cohort workflow.
"""

from os.path import exists
from pathlib import Path

import toml
from cpg_utils import to_path
from pytest_mock import MockFixture

from . import set_config, update_dict

ref_prefix = to_path(__file__).parent / 'data/large_cohort/reference'
gnomad_prefix = ref_prefix / 'gnomad/v0'
broad_prefix = ref_prefix / 'hg38/v0'

TOML = """
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
intervals = ["chr20:start-end", "chrX:start-end", "chrY:start-end"]

[references]
genome_build = "GRCh38"

[storage.default]
default = "{directory}"
web = "{directory}-web"
analysis = "{directory}-analysis"
tmp = "{directory}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.thousand-genomes]
default = "{directory}"
web = "{directory}-web"
analysis = "{directory}-analysis"
tmp = "{directory}-test-tmp"
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

DEFAULT_CONFIG = Path(
    to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml'
)
LARGE_COHORT_CONFIG = Path(
    to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml'
)


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


def test_large_cohort(mocker: MockFixture, tmp_path):
    """
    Run entire workflow in a local mode.
    """
    conf = TOML.format(
        directory=str(tmp_path),
        gnomad_prefix=gnomad_prefix,
        broad_prefix=broad_prefix,
    )
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
    )

    mocker.patch('cpg_workflows.inputs.create_cohort', _mock_cohort)

    from cpg_utils.hail_batch import start_query_context

    from cpg_workflows.large_cohort import (
        ancestry_pca,
        ancestry_plots,
        combiner,
        dense_subset,
        relatedness,
        sample_qc,
        site_only_vcf,
    )

    start_query_context()

    res_pref = tmp_path
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
    ancestry_sample_qc_ht_path = res_pref / 'ancestry_sample_qc.ht'

    ancestry_pca.run(
        dense_mt_path=dense_mt_path,
        sample_qc_ht_path=sample_qc_ht_path,
        relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=res_pref / 'tmp',
        out_scores_ht_path=scores_ht_path,
        out_eigenvalues_ht_path=eigenvalues_ht_path,
        out_loadings_ht_path=loadings_ht_path,
        out_inferred_pop_ht_path=inferred_pop_ht_path,
        out_sample_qc_ht_path=ancestry_sample_qc_ht_path,
    )
    ancestry_plots.run(
        out_path_pattern=res_pref / 'plots' / '{scope}_pc{pci}_{pca_suffix}.{ext}',
        sample_qc_ht_path=ancestry_sample_qc_ht_path,
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
    assert exists(res_pref / 'plots' / 'dataset_pc1_hgdp_1kg_sites.html')
    assert exists(siteonly_vcf_path)
