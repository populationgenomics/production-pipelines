import re
from pathlib import Path
from functools import cached_property

from cpg_workflows.large_cohort import ancestry_pca
from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

from .. import set_config

from .helpers import get_command_str

from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.sequencing_group import create_sequencing_group


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='ancestry-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'dragmap': 'dragmap_image:1.3.0',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                }
            },
        },
        large_cohort={
            'min_pop_prob': 0.5,
            'n_pcs': 16,
            'training_pop': 'Superpopulation name',
        },
    )


class TestAncestryPCA:
    def test_ancestry_pca_typical_run(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test
        batch = create_local_batch(tmp_path)
        # j = dataproc_job(
        #     job_name=self.__class__.__name__,
        #     function=ancestry_pca.run,
        #     function_path_args=dict(
        #         dense_mt_path=tmp_path / 'DenseSubset.mt',
        #         sample_qc_ht_path=tmp_path / 'SampleQC.ht',
        #         relateds_to_drop_ht_path=tmp_path / 'Relatedness.mt',
        #         tmp_prefix=tmp_path,
        #         out_scores_ht_path=tmp_path / 'outputs' / 'scores.ht',
        #         out_eigenvalues_ht_path=tmp_path / 'outputs' / 'eigenvalues.ht',
        #         out_loadings_ht_path=tmp_path / 'outputs' / 'loadings.ht',
        #         out_inferred_pop_ht_path=tmp_path / 'outputs' / 'inferred_pop.ht',
        #         out_sample_qc_ht_path=tmp_path / 'outputs' / 'sample_qc.ht',
        #     ),
        #     # depends_on=inputs.get_jobs(cohort),
        # )

        # ---- Assertions
        assert False

    def test_add_background(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test

        # ---- Assertions

        assert False

    def test_fail_n_pcs_less_than_min_n_pcs(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test

        # ---- Assertions

        assert False

    def test_run_pca_ancestry_analysis(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test

        # ---- Assertions

        assert False

    def test_infer_pop_labels(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test

        # ---- Assertions

        assert False

    def test_run_assign_population_pcs(
        self, tmp_path: Path, pop_pca_scores_ht_, min_prob_
    ):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')

        # ---- The job that we want to test

        # ---- Assertions

        assert False
