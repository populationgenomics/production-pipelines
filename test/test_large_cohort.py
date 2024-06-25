"""
Test large-cohort workflow.
"""

import os
from os.path import exists
from pathlib import Path

from pytest_mock import MockFixture

from cpg_utils import Path as CPGPath
from cpg_utils import to_path
from cpg_utils.hail_batch import start_query_context
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.large_cohort import (
    ancestry_pca,
    ancestry_plots,
    combiner,
    dense_subset,
    relatedness,
    sample_qc,
    site_only_vcf,
)
from cpg_workflows.targets import Cohort

from . import set_config
from .factories.config import HailConfig, PipelineConfig, StorageConfig, WorkflowConfig
from .factories.types import SequencingType

ref_prefix = to_path(__file__).parent / 'data/large_cohort/reference'
gnomad_prefix = ref_prefix / 'gnomad/v0'
broad_prefix = ref_prefix / 'hg38/v0'


DEFAULT_CONFIG = Path(to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml')
LARGE_COHORT_CONFIG = Path(to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml')


def create_config(
    tmp_path: Path,
    seq_type: SequencingType,
    gnomad_prefix: CPGPath,
    broad_prefix: CPGPath,
) -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='large-cohort-test',
            access_level='test',
            sequencing_type=seq_type,
            check_intermediates=True,
        ),
        hail=HailConfig(query_backend='spark_local'),
        storage={
            'default': StorageConfig(
                default=tmp_path,
                web=tmp_path / 'web',
                analysis=tmp_path / 'analysis',
                tmp=tmp_path / 'test-tmp',
            ),
            'large-cohort-test': StorageConfig(
                default=tmp_path,
                web=tmp_path / 'web',
                analysis=tmp_path / 'analysis',
                tmp=tmp_path / 'test-tmp',
            ),
        },
        references={
            'genome_build': 'GRCh38',
            'ancestry': {
                'sites_table': (gnomad_prefix / 'sample_qc-test' / 'pre_ld_pruning_qc_variants.ht'),
            },
            'gnomad': {
                'tel_and_cent_ht': (
                    gnomad_prefix / 'telomeres_and_centromeres' / 'hg38.telomeresAndMergedCentromeres.ht'
                ),
                'lcr_intervals_ht': (gnomad_prefix / 'lcr_intervals' / 'LCRFromHengHg38.ht'),
                'seg_dup_intervals_ht': (gnomad_prefix / 'seg_dup_intervals' / 'GRCh38_segdups.ht'),
                'clinvar_ht': (gnomad_prefix / 'clinvar' / 'clinvar_20190923.ht'),
                'hapmap_ht': (gnomad_prefix / 'hapmap' / 'hapmap_3.3.hg38.ht'),
                'kgp_omni_ht': (gnomad_prefix / 'kgp' / '1000G_omni2.5.hg38.ht'),
                'kgp_hc_ht': (gnomad_prefix / 'kgp' / '1000G_phase1.snps.high_confidence.hg38.ht'),
                'mills_ht': (gnomad_prefix / 'mills' / 'Mills_and_1000G_gold_standard.indels.hg38.ht'),
            },
            'gatk_sv': {
                'protein_coding_gtf': (
                    broad_prefix / 'sv-resources' / 'resources' / 'v1' / 'MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
                ),
            },
            'broad': {
                'genome_calling_interval_lists': (broad_prefix / 'wgs_calling_regions.hg38.interval_list'),
            },
        },
        large_cohort={
            'n_pcs': 3,
            'sample_qc_cutoffs': {
                'min_n_snps': 2500,
            },
            'combiner': {'intervals': ['chr20:start-end', 'chrX:start-end', 'chrY:start-end']},
        },
    )


def _mock_cohort(dataset_id: str):
    cohort = Cohort()
    dataset = cohort.create_dataset(dataset_id)
    gvcf_root = to_path(__file__).parent / 'data' / 'large_cohort' / 'gvcf'
    found_gvcf_paths = list(gvcf_root.glob('*.g.vcf.gz'))
    assert len(found_gvcf_paths) > 0, gvcf_root
    for gvcf_path in found_gvcf_paths:
        sequencing_group_id = gvcf_path.name.split('.')[0]
        sequencing_group = dataset.add_sequencing_group(
            id=sequencing_group_id,
            external_id=sequencing_group_id.replace('CPG', 'EXT'),
        )
        sequencing_group.gvcf = GvcfPath(gvcf_path)
    return cohort


class TestAllLargeCohortMethods:
    def test_with_sample_data(self, mocker: MockFixture, tmp_path: Path):
        """
        Run entire workflow in a local mode.
        """
        conf = create_config(tmp_path, 'genome', gnomad_prefix, broad_prefix)
        set_config(
            conf,
            tmp_path / 'config.toml',
            merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
        )

        mocker.patch(
            'cpg_workflows.inputs.deprecated_create_cohort',
            lambda: _mock_cohort(conf.workflow.dataset),
        )
        # skip can_reuse, implicit skip of existence checks
        mocker.patch('cpg_workflows.large_cohort.combiner.can_reuse', lambda x: False)
        mocker.patch('cpg_workflows.large_cohort.relatedness.can_reuse', lambda x: False)
        mocker.patch('cpg_workflows.large_cohort.site_only_vcf.can_reuse', lambda x: False)

        start_query_context()

        res_pref = tmp_path
        vds_path = res_pref / 'v01.vds'
        combiner.run(out_vds_path=vds_path, tmp_prefix=res_pref / 'tmp')

        sample_qc_ht_path = res_pref / 'sample_qc.ht'
        sample_qc.run(
            vds_path=str(vds_path),
            out_sample_qc_ht_path=str(sample_qc_ht_path),
            tmp_prefix=os.path.join(res_pref, 'tmp'),
        )

        dense_mt_path = res_pref / 'dense.mt'
        dense_subset.run(
            vds_path=str(vds_path),
            out_dense_mt_path=str(dense_mt_path),
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
            relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        )

        siteonly_vcf_path = res_pref / 'siteonly.vcf.bgz'
        site_only_vcf.run(
            vds_path=str(vds_path),
            sample_qc_ht_path=str(sample_qc_ht_path),
            relateds_to_drop_ht_path=str(relateds_to_drop_ht_path),
            out_vcf_path=str(siteonly_vcf_path),
            tmp_prefix=str(res_pref / 'tmp'),
        )

        assert exists(vds_path)
        assert exists(res_pref / 'plots' / 'dataset_pc1_hgdp_1kg_sites.html')
        assert exists(siteonly_vcf_path)
