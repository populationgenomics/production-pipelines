#!/usr/bin/env python3

"""
Stages to run somalier tools.
"""

import logging

from cpg_pipes.jobs import somalier
from cpg_pipes.pipeline.analysis import CramPath
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import DatasetStage, StageInput, StageOutput, SampleStage

logger = logging.getLogger(__file__)


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
        cram_path = sample.get_cram_path()

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


@stage(required_stages=CramSomalierStage)
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
        fp_by_sid = {
            s.id: inputs.as_path(stage=CramSomalierStage, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_result(dataset)['html']
        html_url = str(html_path).replace(
            str(dataset.get_web_bucket()), dataset.get_web_url()
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


@stage(required_stages=CramSomalierStage)
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
        fp_by_sid = {
            s.id: inputs.as_path(stage=CramSomalierStage, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_result(dataset)['html']
        html_url = str(html_path).replace(
            str(dataset.get_web_bucket()), dataset.get_web_url()
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
