"""
Test large-cohort workflow.
"""

from os.path import join as os_path_join
from pathlib import Path
from zipfile import ZipFile

from pytest_mock import MockFixture

import cpg_workflows.inputs
from cpg_utils import Path as CPGPath
from cpg_utils import to_path
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
from cpg_workflows.targets import MultiCohort
from cpg_workflows.utils import exists_not_cached

from . import set_config
from .factories.config import HailConfig, PipelineConfig, StorageConfig, WorkflowConfig
from .factories.types import SequencingType

ref_prefix = to_path(__file__).parent / 'data/large_cohort/reference'

# input files for each component
compressed_inputs = to_path(__file__).parent / 'data/large_cohort/compressed_dirs'
compressed_vds_path = str(compressed_inputs / 'v01.vds.zip')
compressed_dense_mt_path = str(compressed_inputs / 'dense.mt.zip')
compressed_sample_qc_ht_path = str(compressed_inputs / 'sample_qc.ht.zip')
compressed_relateds_to_drop_ht_path = str(compressed_inputs / 'relateds_to_drop.ht.zip')
# compressed_relatedness_ht_path = str(compressed_inputs / 'relatedness.ht.zip')
compressed_loadings_ht_path = str(compressed_inputs / 'loadings.ht.zip')
compressed_scores_ht_path = str(compressed_inputs / 'scores.ht.zip')
compressed_eigen_ht_path = str(compressed_inputs / 'eigenvalues.ht.zip')
compressed_inferred_pop_ht_path = str(compressed_inputs / 'inferred_pop.ht.zip')
compressed_ancestry_sample_qc_ht_path = str(compressed_inputs / 'ancestry_sample_qc.ht.zip')

GNOMAD_PREFIX = ref_prefix / 'gnomad/v0'
BROAD_PREFIX = ref_prefix / 'hg38/v0'

DEFAULT_CONFIG = Path(to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml')
LARGE_COHORT_CONFIG = Path(to_path(__file__).parent.parent / 'configs' / 'defaults' / 'large_cohort.toml')


def create_config(
    tmp_path: Path,
    sequencing_type: SequencingType = 'genome',
    gnomad_prefix: CPGPath = GNOMAD_PREFIX,
    broad_prefix: CPGPath = BROAD_PREFIX,
) -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='large-cohort-test',
            access_level='test',
            sequencing_type=sequencing_type,
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
            'combiner': {
                'intervals': ['chr20:start-end', 'chrX:start-end', 'chrY:start-end'],
                'force_new_combiner': False,
            },
        },
    )


def _mock_cohort(dataset_id: str):
    mc = MultiCohort()
    cohort = mc.create_cohort(id='COH1', name='large-cohort-test')
    dataset = cohort.create_dataset(dataset_id)
    mc_dataset = mc.add_dataset(dataset)
    gvcf_root = to_path(__file__).parent / 'data' / 'large_cohort' / 'gvcf'
    found_gvcf_paths = list(gvcf_root.glob('*.g.vcf.gz'))
    assert len(found_gvcf_paths) > 0, gvcf_root
    for gvcf_path in found_gvcf_paths:
        sequencing_group_id = gvcf_path.name.split('.')[0]
        sequencing_group = dataset.add_sequencing_group(
            id=sequencing_group_id,
            external_id=sequencing_group_id.replace('CPG', 'EXT'),
            sequencing_type='genome',
            sequencing_technology='short-read',
            sequencing_platform='illumina',
        )
        sequencing_group.gvcf = GvcfPath(gvcf_path)
        mc_dataset.add_sequencing_group_object(sequencing_group)
    return mc


def decompress_into_job_tmp(tmp_path: Path, compressed_paths: list[str]):
    """
    Takes a list of compressed paths, and decompresses them into a job temp location
    """
    for path in compressed_paths:
        with ZipFile(path, 'r') as zip_ref:
            zip_ref.extractall(tmp_path)


def test_combiner(mocker: MockFixture, tmp_path: Path):
    conf = create_config(tmp_path)
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

    test_cohort: MultiCohort = cpg_workflows.inputs.deprecated_create_cohort()
    gvcf_paths: list[str] = [str(sg.gvcf) for sg in test_cohort.get_sequencing_groups()]

    vds_path = str(tmp_path / 'v01.vds')

    # we're passing a specific minority of intervals here, to test that the combiner works on a timely test case
    combiner.run(
        output_vds_path=vds_path,
        sequencing_type=conf['workflow']['sequencing_type'],
        tmp_prefix=str(tmp_path / 'tmp'),
        genome_build=conf['references']['genome_build'],
        gvcf_paths=gvcf_paths,
        vds_paths=None,
        save_path=None,
        specific_intervals=conf['large_cohort']['combiner']['intervals'],
        force_new_combiner=conf['large_cohort']['combiner']['force_new_combiner'],
    )

    # do some testing here
    assert exists_not_cached(vds_path)


def test_sample_qc(mocker: MockFixture, tmp_path: Path):
    conf = create_config(tmp_path)
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
    mocker.patch('cpg_workflows.large_cohort.sample_qc.can_reuse', lambda x: False)

    # open that VDS into a job temp location
    decompress_into_job_tmp(tmp_path, [compressed_vds_path])

    sample_qc_ht_path = tmp_path / 'sample_qc.ht'
    sample_qc.run(
        vds_path=str(os_path_join(tmp_path, 'v01.vds')),
        out_sample_qc_ht_path=str(sample_qc_ht_path),
        tmp_prefix=os_path_join(tmp_path, 'tmp'),
    )

    # do some testing here
    assert exists_not_cached(sample_qc_ht_path)


def test_densify_mt(tmp_path: Path):
    conf = create_config(tmp_path)
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
    )

    # open that VDS into a job temp location
    decompress_into_job_tmp(tmp_path, [compressed_vds_path])

    dense_mt_output_path = tmp_path / 'dense.mt'

    # uses get_config()['references']['ancestry']['sites_table']
    dense_subset.run(
        vds_path=str(tmp_path / 'v01.vds'),
        out_dense_mt_path=str(dense_mt_output_path),
    )

    # do some testing here
    assert exists_not_cached(dense_mt_output_path)


def test_relatedness(mocker: MockFixture, tmp_path: Path):
    conf = create_config(tmp_path)
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
    )

    # skip can_reuse, implicit skip of existence checks
    mocker.patch('cpg_workflows.large_cohort.relatedness.can_reuse', lambda x: False)

    # decompress the sample QC HT and dense MT
    decompress_into_job_tmp(tmp_path, [compressed_sample_qc_ht_path, compressed_dense_mt_path])

    relateds_to_drop_ht_path = tmp_path / 'relateds_to_drop.ht'
    relatedness_ht_path = tmp_path / 'relatedness.ht'

    # uses get_config()['large_cohort']['remove_failed_qc_pca'] and get_config()['large_cohort']['max_kin']
    relatedness.run(
        dense_mt_path=tmp_path / 'dense.mt',
        sample_qc_ht_path=tmp_path / 'sample_qc.ht',
        out_relatedness_ht_path=relatedness_ht_path,
        out_relateds_to_drop_ht_path=relateds_to_drop_ht_path,
        tmp_prefix=tmp_path / 'tmp',
    )

    # do some testing here
    assert exists_not_cached(relateds_to_drop_ht_path)
    assert exists_not_cached(relatedness_ht_path)


def test_site_only(mocker: MockFixture, tmp_path: Path):
    # skip can_reuse, implicit skip of existence checks
    mocker.patch('cpg_workflows.large_cohort.site_only_vcf.can_reuse', lambda x: False)

    # decompress the sample QC HT, VDS, and relateds to drop HT
    decompress_into_job_tmp(
        tmp_path,
        [compressed_sample_qc_ht_path, compressed_vds_path, compressed_relateds_to_drop_ht_path],
    )

    siteonly_vcf_path = tmp_path / 'siteonly.vcf.bgz'
    site_only_vcf.run(
        vds_path=str(tmp_path / 'v01.vds'),
        sample_qc_ht_path=str(tmp_path / 'sample_qc.ht'),
        relateds_to_drop_ht_path=str(tmp_path / 'relateds_to_drop.ht'),
        out_vcf_path=str(siteonly_vcf_path),
        tmp_prefix=str(tmp_path / 'tmp'),
    )

    # do some testing here
    assert exists_not_cached(siteonly_vcf_path)


def test_ancestry(tmp_path: Path):
    conf = create_config(tmp_path)
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
    )

    scores_ht_path = tmp_path / 'scores.ht'
    eigenvalues_ht_path = tmp_path / 'eigenvalues.ht'
    loadings_ht_path = tmp_path / 'loadings.ht'
    inferred_pop_ht_path = tmp_path / 'inferred_pop.ht'
    ancestry_sample_qc_ht_path = tmp_path / 'ancestry_sample_qc.ht'

    # decompress the dense MT, sample_qc, and relateds to drop HT
    decompress_into_job_tmp(
        tmp_path,
        [
            compressed_dense_mt_path,
            compressed_sample_qc_ht_path,
            compressed_relateds_to_drop_ht_path,
        ],
    )

    # uses get_config()['large_cohort']['min_pop_prob'] & get_config()['large_cohort'].get('pca_samples_to_remove', [])
    ancestry_pca.run(
        dense_mt_path=tmp_path / 'dense.mt',
        sample_qc_ht_path=tmp_path / 'sample_qc.ht',
        relateds_to_drop_ht_path=tmp_path / 'relateds_to_drop.ht',
        tmp_prefix=tmp_path / 'tmp',
        out_scores_ht_path=scores_ht_path,
        out_eigenvalues_ht_path=eigenvalues_ht_path,
        out_loadings_ht_path=loadings_ht_path,
        out_inferred_pop_ht_path=inferred_pop_ht_path,
        out_sample_qc_ht_path=ancestry_sample_qc_ht_path,
    )

    # do some testing here
    for output in [
        scores_ht_path,
        eigenvalues_ht_path,
        loadings_ht_path,
        inferred_pop_ht_path,
        ancestry_sample_qc_ht_path,
    ]:
        assert exists_not_cached(output)


def test_ancestry_plots(tmp_path: Path):
    conf = create_config(tmp_path)
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, LARGE_COHORT_CONFIG],
    )

    # decompress all the inputs tables for plotting ancestry
    decompress_into_job_tmp(
        tmp_path,
        [
            compressed_ancestry_sample_qc_ht_path,
            compressed_scores_ht_path,
            compressed_eigen_ht_path,
            compressed_loadings_ht_path,
            compressed_inferred_pop_ht_path,
            compressed_relateds_to_drop_ht_path,
        ],
    )

    # uses a few config entries
    ancestry_plots.run(
        out_path_pattern=tmp_path / 'plots' / '{scope}_pc{pci}_{pca_suffix}.{ext}',
        sample_qc_ht_path=tmp_path / 'ancestry_sample_qc.ht',
        scores_ht_path=tmp_path / 'scores.ht',
        eigenvalues_ht_path=tmp_path / 'eigenvalues.ht',
        loadings_ht_path=tmp_path / 'loadings.ht',
        inferred_pop_ht_path=tmp_path / 'inferred_pop.ht',
        relateds_to_drop_ht_path=tmp_path / 'relateds_to_drop.ht',
    )

    # do some testing here
    assert exists_not_cached(tmp_path / 'plots' / 'dataset_pc1_hgdp_1kg_sites.html')
