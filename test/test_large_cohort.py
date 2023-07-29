"""
Test large-cohort workflow.
"""

from os.path import exists
from pathlib import Path

import pytest
from cpg_utils import Path as CPGPath
from cpg_utils import to_path
from cpg_utils.hail_batch import start_query_context
from pytest_mock import MockFixture

from hail.utils import FatalError

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
from cpg_workflows.targets import Cohort, Dataset

from . import set_config
from .factories.config import HailConfig, PipelineConfig, StorageConfig, WorkflowConfig
from .factories.sequencing_group import create_sequencing_group
from .factories.types import SequencingType
# test imports
import pytest
from hail.utils.java import FatalError

ref_prefix = to_path(__file__).parent / 'data/large_cohort/reference'
gnomad_prefix = ref_prefix / 'gnomad/v0'
broad_prefix = ref_prefix / 'hg38/v0'


DEFAULT_CONFIG = Path(
    to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml'
)
LARGE_COHORT_CONFIG = Path(
    to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml'
)


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
            'gnomad': {
                'tel_and_cent_ht': (
                    gnomad_prefix
                    / 'telomeres_and_centromeres'
                    / 'hg38.telomeresAndMergedCentromeres.ht'
                ),
                'lcr_intervals_ht': (
                    gnomad_prefix / 'lcr_intervals' / 'LCRFromHengHg38.ht'
                ),
                'seg_dup_intervals_ht': (
                    gnomad_prefix / 'seg_dup_intervals' / 'GRCh38_segdups.ht'
                ),
                'clinvar_ht': (gnomad_prefix / 'clinvar' / 'clinvar_20190923.ht'),
                'hapmap_ht': (gnomad_prefix / 'hapmap' / 'hapmap_3.3.hg38.ht'),
                'kgp_omni_ht': (gnomad_prefix / 'kgp' / '1000G_omni2.5.hg38.ht'),
                'kgp_hc_ht': (
                    gnomad_prefix / 'kgp' / '1000G_phase1.snps.high_confidence.hg38.ht'
                ),
                'mills_ht': (
                    gnomad_prefix
                    / 'mills'
                    / 'Mills_and_1000G_gold_standard.indels.hg38.ht'
                ),
                'predetermined_qc_variants': (
                    gnomad_prefix / 'sample_qc' / 'pre_ld_pruning_qc_variants.ht'
                ),
            },
            'broad': {
                'genome_calling_interval_lists': (
                    broad_prefix / 'wgs_calling_regions.hg38.interval_list'
                ),
                'protein_coding_gtf': (
                    broad_prefix
                    / 'sv-resources'
                    / 'resources'
                    / 'v1'
                    / 'MANE.GRCh38.v0.95.select_ensembl_genomic.gtf'
                ),
            },
        },
        large_cohort={
            'n_pcs': 3,
            'sample_qc_cutoffs': {
                'min_n_snps': 2500,
            },
            'combiner': {
                'intervals': ['chr20:start-end', 'chrX:start-end', 'chrY:start-end']
            },
        },
    )


def _mock_cohort(dataset_id: str):
    dataset = Dataset(name=dataset_id)

    # Parse gVCF files from the test/data/large_cohort/gvcf directory
    gvcf_root = to_path(__file__).parent / 'data' / 'large_cohort' / 'gvcf'
    found_gvcf_paths = list(gvcf_root.glob('*.g.vcf.gz'))
    assert len(found_gvcf_paths) > 0, gvcf_root

    for gvcf_path in found_gvcf_paths:
        sequencing_group_id = gvcf_path.name.split('.')[0]
        create_sequencing_group(
            id=sequencing_group_id,
            external_id=sequencing_group_id.replace('CPG', 'EXT'),
            dataset=dataset,
            gvcf=GvcfPath(gvcf_path),
        )

    cohort = Cohort()
    cohort.add_dataset(dataset)
    return cohort

class TestCombiner:
    def test_fails_if_given_invalid_chromosome_that_does_not_exist(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # Creating a config and modyfing it in some way
        conf = create_config(
            tmp_path,
            seq_type='genome',
            gnomad_prefix=gnomad_prefix,
            broad_prefix=broad_prefix,
        )
        conf.large_cohort['combiner']['intervals'] = ['chr27:start-end']

        # Set config and patch cohort
        set_config(
            conf,
            tmp_path / 'config.toml',
            merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
        )

        mocker.patch(
            'cpg_workflows.inputs.create_cohort',
            lambda: _mock_cohort(conf.workflow.dataset),
        )

        # Run the combiner function
        start_query_context()
        res_pref = tmp_path
        vds_path = res_pref / 'v01.vds'
        with pytest.raises(FatalError, match='invalid interval expression'):
            combiner.run(out_vds_path=vds_path, tmp_prefix=res_pref / 'tmp')

    @MockFixture
    def test_uses_default_genome_intervals_if_intervals_are_not_specified(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # Test that if conf.large_cohort['combiner']['intervals'] is empty
        # then hl.vds.new_combiner is called with use_genome_default_intervals=True

        # Creating a config and modyfing it in some way
        conf = create_config(
            tmp_path,
            seq_type='genome',
            gnomad_prefix=gnomad_prefix,
            broad_prefix=broad_prefix,
        )
        conf.large_cohort['combiner']['intervals'] = [] # set to empty list

        # Set config and patch cohort
        set_config(
            conf,
            tmp_path / 'config.toml',
            merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
        )

        mocker.patch(
            'cpg_workflows.inputs.create_cohort',
            lambda: _mock_cohort(conf.workflow.dataset),
        )

        # Run the combiner function
        start_query_context()
        res_pref = tmp_path
        vds_path = res_pref / 'v01.vds'
        combined = combiner.run(out_vds_path=vds_path, tmp_prefix=res_pref / 'tmp')
        print()
        # wrtie a pytest_mock that replicates a combiner.run() call and tests whether
        # 'params' attribute = {'intervals': [], 'use_genome_default_intervals': True}
        # Hint: read about pytest_mock spy from the guide in the testing google document
        

    def test_uses_default_exome_intervals_if_intervals_are_not_specified(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # Test that if conf.large_cohort['combiner']['intervals'] is empty
        # then hl.vds.new_combiner is called with use_exome_default_intervals=True

        # Hint: read about pytest_mock spy from the guide in the testing google document
        pass

    def test_fails_if_given_malformed_intervals(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # See https://hail.is/docs/0.2/functions/genetics.html#hail.expr.functions.parse_locus_interval
        # for valid formats
        pass

    def test_fails_if_given_duplicate_sequencing_groups(
        self, mocker: MockFixture, tmp_path: Path
    ):
        pass

    def test_fails_if_given_two_sequencing_groups_with_the_same_gvcf_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        pass

    @pytest.mark.parametrize('seq_type', ['exome', 'genome'])
    def test_calls_hail_combiner_with_correct_parameters(
        self, mocker: MockFixture, tmp_path: Path, seq_type: SequencingType
    ):
        pass

    def test_fails_if_all_sequencing_groups_do_not_have_a_gvcf_file(
        self, mocker: MockFixture, tmp_path: Path
    ):
        pass

    def test_can_reuse_existing_vds_that_exists_at_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # Hint: create a blank file and mock hl.vds.read_vds to fake a return value
        pass


class TestAllLargeCohortMethods:
    def test_with_sample_data(self, mocker: MockFixture, tmp_path: Path):
        """
        Run entire workflow in a local mode.
        """
        conf = create_config(
            tmp_path,
            seq_type='genome',
            gnomad_prefix=gnomad_prefix,
            broad_prefix=broad_prefix,
        )

        set_config(
            conf,
            tmp_path / 'config.toml',
            merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
        )

        mocker.patch(
            'cpg_workflows.inputs.create_cohort',
            lambda: _mock_cohort(conf.workflow.dataset),
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
