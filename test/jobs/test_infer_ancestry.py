import re
from pathlib import Path
from functools import cached_property

from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
from cpg_workflows.large_cohort.ancestry_pca import add_background, run, _run_pca_ancestry_analysis, _infer_pop_labels

from .. import set_config

from .helpers import get_command_str

from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.sequencing_group import create_sequencing_group

class TestAncestryPCA:
    @cached_property
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
                'training_pop': 'Superpopulation name'
            }
        )
    
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
    
    def test_typical_run(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')
        
        # ---- The job that we want to test
        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort, Relatedness, key='relateds_to_drop'
                ),
                tmp_prefix=self.tmp_prefix,
                out_scores_ht_path=self.expected_outputs(cohort)['scores'],
                out_eigenvalues_ht_path=self.expected_outputs(cohort)['eigenvalues'],
                out_loadings_ht_path=self.expected_outputs(cohort)['loadings'],
                out_inferred_pop_ht_path=self.expected_outputs(cohort)['inferred_pop'],
                out_sample_qc_ht_path=self.expected_outputs(cohort)['sample_qc_ht'],
            ),
            depends_on=inputs.get_jobs(cohort),
        )

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
    
    def test_run_assign_population_pcs(self, tmp_path: Path, pop_pca_scores_ht_, min_prob_):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')
        
        # ---- The job that we want to test

        # ---- Assertions

        assert False
    

    def test_creates_one_align_job(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config, tmp_path / 'config.toml')
        


        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)

        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
        )

        # ---- The job that we want to test
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        assert len(align_jobs) == 1

    def test_sorts_output_with_bamtools(self, tmp_path: Path):
        # ---- Test setup
        config = self.default_config()
        set_config(config.as_dict(), tmp_path / 'config.toml')
        


        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)

        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
        )

        # ---- The job that we want to test
        _ = align(
            b=batch,
            sequencing_group=sg,
            extra_label=dataset_id,
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.NO_MARKDUP,
        )

        # ---- Assertions
        align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
        cmd = get_command_str(align_jobs[0])
        assert re.search(r'\| samtools sort .* -Obam', cmd)
