"""
Stage to make a TrainGCNV model. We don't want to make it on each workflow run.
"""
import os

import click

from cpg_utils import Path
from cpg_utils.config import set_config_paths
from cpg_utils.workflows.targets import Cohort
from cpg_utils.workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    run_workflow,
)
from .gatk_sv import GatherSampleEvidence
from .gatk_sv import add_gatk_sv_job, get_references, get_images, SV_CALLERS


@stage(required_stages=[GatherSampleEvidence])
class TrainGCNV(CohortStage):
    """
    # https://github.com/populationgenomics/gatk-sv/blob/main/wdl/TrainGCNV.wdl
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to return 2 tarballs
        """
        d = dict()
        fname_by_key = {
            'cohort_contig_ploidy_model_tar': 'cohort_contig_ploidy_model.tar',
            'cohort_contig_ploidy_calls_tar': 'cohort_contig_ploidy_calls.tar',
        }
        for key, fname in fname_by_key.items():
            d[key] = (
                cohort.analysis_dataset.get_bucket()
                / 'gatk_sv'
                / self.name.lower()
                / fname
            )
        return d

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        d = inputs.as_dict_by_target(GatherSampleEvidence)

        sids = cohort.get_sample_ids()
        input_dict = {
            'cohort': cohort.target_id,
            'samples': sids,
            'count_files': [str(d[sid]['coverage_counts']) for sid in sids],
            'ref_copy_number_autosomal_contigs': 2,
            'num_intervals_per_scatter': 5000,
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(d[sid][f'{caller}_vcf']) for sid in sids
            ]

        input_dict.update(
            get_images(
                [
                    'sv_base_mini_docker',
                    'condense_counts_docker',
                    'gatk_docker',
                    'linux_docker',
                ]
            )
        )

        input_dict.update(
            get_references(
                [
                    'reference_fasta',
                    'reference_index',
                    'reference_dict',
                    'allosomal_contigs',
                    'contig_ploidy_priors',
                ]
            )
        )

        expected_d = self.expected_outputs(cohort)

        j = add_gatk_sv_job(
            batch=self.b,
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )

        return self.make_outputs(cohort, data=expected_d, jobs=[j])


@click.command()
@click.argument('config_paths', nargs=-1)
def main(config_paths: list[str]):
    """
    Run a workflow, using CONFIG_PATHS in the order specified, overriding
    $CPG_CONFIG_PATH if specified.
    """
    if _cpg_config_path_env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths = _cpg_config_path_env_var.split(',') + config_paths
    set_config_paths(list(config_paths))
    run_workflow(stages=[TrainGCNV])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
