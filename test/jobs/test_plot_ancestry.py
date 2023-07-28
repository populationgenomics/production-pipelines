import re
import git
from pathlib import Path
from functools import cached_property

from cpg_workflows.large_cohort import ancestry_plots

from .. import set_config

from .helpers import get_command_str

from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.sequencing_group import create_sequencing_group


class TestAncestryPlots:
    @cached_property
    def default_config(self) -> PipelineConfig:
        return PipelineConfig(
            workflow=WorkflowConfig(
                dataset='ancestry-test',
                access_level='test',
                sequencing_type='genome',
                check_inputs=False,
                scatter_count=20,
                dataset_gcp_project='cpg-fake-gcp-project',
                driver_image='fake-driver-image',
            ),
            images={
                'dragmap': 'dragmap_image:1.3.0',
            },
            references={
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                }
            },
            large_cohort={
                'min_pop_prob': 0.5,
                'n_pcs': 16,
                'training_pop': 'Superpopulation name',
            },
        )

    def test_ancestry_plots_typical_run(self, tmp_path: Path):
        # ---- Test setup
        set_config(self.default_config, tmp_path / 'config.toml')

        # TODO: move to correct location after analysis_runner bug fixed
        # see: https://github.com/orgs/populationgenomics/projects/17/views/1?pane=issue&itemId=34321186
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        # ---- The job that we want to test
        # -- args
        out_path_pattern = tmp_path
        sample_qc_ht_path = tmp_path / 'SampleQC.ht'
        scores_ht_path = tmp_path / 'Scores.ht'
        eigenvalues_ht_path = tmp_path / 'eigenvalues.ht'
        loadings_ht_path = tmp_path / 'loading.ht'
        inferred_pop_ht_path = tmp_path / 'inferred_pop.ht'

        # -- create job
        _ = create_local_batch(tmp_path)
        job = dataproc_job(
            job_name=self.__class__.__name__,
            function=ancestry_plots.run,
            function_path_args=dict(
                out_path_pattern=out_path_pattern,
                sample_qc_ht_path=sample_qc_ht_path,
                scores_ht_path=scores_ht_path,
                eigenvalues_ht_path=eigenvalues_ht_path,
                loadings_ht_path=loadings_ht_path,
                inferred_pop_ht_path=inferred_pop_ht_path,
            ),
        )
        cmd = get_command_str(job)
        repo = git.Repo(search_parent_directories=True)
        git_hash = repo.head.object.hexsha

        # ---- Assertions
        assert re.search(git_hash, cmd)
        assert re.search('dataproc *submit *--region=australia-southeast1', cmd)
        assert re.search('--pyfiles *cpg_workflows,gnomad_methods/gnomad', cmd)
        assert re.search('cpg_workflows/large_cohort/dataproc_script.py', cmd)
        assert re.search('cpg_workflows.large_cohort.ancestry_plots', cmd)
        assert re.search(f'-p.*{out_path_pattern}', cmd)
        assert re.search(f'-p.*{sample_qc_ht_path}', cmd)
        assert re.search(f'-p.*{scores_ht_path}', cmd)
        assert re.search(f'-p.*{eigenvalues_ht_path}', cmd)
        assert re.search(f'-p.*{loadings_ht_path}', cmd)
        assert re.search(f'-p.*{inferred_pop_ht_path}', cmd)
