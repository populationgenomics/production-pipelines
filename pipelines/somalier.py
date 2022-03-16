#!/usr/bin/env python3

"""
Batch pipeline to check pedigree on samples.
"""

import logging

import click

from cpg_pipes.jobs import somalier
from cpg_pipes.pipeline.analysis import CramPath
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import DatasetStage, StageInput, StageOutput, SampleStage
from cpg_pipes.pipeline.pipeline import stage, Pipeline, PipelineError
from cpg_pipes.pipeline.cli_opts import pipeline_click_options

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class CramSomalierStage(SampleStage):
    """
    Genereate fingerprints from CRAMs for pedigree checks.
    """

    def expected_result(self, sample: Sample):
        """
        Expected to generate the fingerprints file
        """
        return sample.get_cram_path().somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using a function from the jobs module.
        """
        cram_path = sample.analysis_cram_path()
        if not cram_path:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(
                    f'Could not find CRAM analysis, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(f'No CRAM analysis found for {sample.id}')

        expected_path = self.expected_result(sample)
        j = somalier.extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=CramPath(cram_path),
            out_fpath=expected_path,
            overwrite=not self.pipe.check_intermediates,
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage(required_stages=CramSomalierStage, forced=True)
class CramSomalierPedigree(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints.
    """

    def expected_result(self, dataset: Dataset):
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """

        prefix = dataset.get_analysis_bucket() / 'qc' / 'somalier'
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': dataset.get_web_bucket() / 'pedigree' / 'somalier.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = inputs.as_path_by_target(stage=CramSomalierStage)

        html_path = self.expected_result(dataset)['html']
        html_url = str(html_path).replace(
            str(html_path.parent), dataset.get_web_url()
        )

        j = somalier.pedigree(
            self.pipe.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.pipe.check_intermediates,
            out_samples_path=self.expected_result(dataset)['samples'],
            out_pairs_path=self.expected_result(dataset)['pairs'],
            out_html_path=html_path,
            out_html_url=html_url,
            depends_on=inputs.get_jobs(),
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


@stage(required_stages=CramSomalierStage, forced=True)
class CramSomalierAncestry(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints
    """

    def expected_result(self, dataset: Dataset):
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        *.somalier-ancestry.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """

        prefix = dataset.get_analysis_bucket() / 'ancestry'
        return {
            'tsv': prefix / f'{dataset.name}.somalier-ancestry.tsv',
            'html': dataset.get_web_bucket() / 'somalier-ancestry.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = inputs.as_path_by_target(stage=CramSomalierStage)
        
        html_path = self.expected_result(dataset)['html']
        html_url = str(html_path).replace(
            str(html_path.parent), dataset.get_web_url()
        )

        j = somalier.ancestry(
            self.pipe.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.pipe.check_intermediates,
            out_tsv_path=self.expected_result(dataset)['tsv'],
            out_html_path=html_path,
            out_html_url=html_url,
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):  # pylint: disable=missing-function-docstring
    Pipeline(
        name='pedigree_check',
        description='Pedigree checks',
        **kwargs,
    ).submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
