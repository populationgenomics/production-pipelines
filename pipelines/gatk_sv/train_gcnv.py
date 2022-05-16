"""
Stage to make a TrainGCNV model. We don't want to make it on each pipeline run,
and mostly a pre-trained 1kg reference is sufficient.
"""

import logging

import click

from cpg_pipes import Path
from cpg_pipes.pipeline import stage, StageOutput, CohortStage, \
    pipeline_click_options, \
    create_pipeline
from cpg_pipes.pipeline.pipeline import StageInput
from cpg_pipes.targets import Cohort

from pipelines.gatk_sv.gatk_sv import GatherSampleEvidence
from pipelines.gatk_sv.utils import add_gatksv_job, get_references, get_dockers, \
    SV_CALLERS

logger = logging.getLogger(__file__)


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
            d[
                key] = cohort.analysis_dataset.get_bucket() / 'gatk_sv' / \
                       self.name.lower() / fname
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
            get_dockers(
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

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset=cohort.analysis_dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )

        return self.make_outputs(cohort, data=output_dict, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):
    """
    GATK-SV workflow.
    """
    pipeline = create_pipeline(
        name='gatk_sv_train_gcnv',
        **kwargs
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
