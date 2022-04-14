#!/usr/bin/env python3

"""
Stages to run somalier tools.
"""

import logging

from .. import Path
from ..targets import Dataset, Sample
from ..pipeline import stage, SampleStage, DatasetStage, StageInput, StageOutput
from ..jobs import somalier

logger = logging.getLogger(__file__)


@stage
class CramSomalier(SampleStage):
    """
    Genereate fingerprints from CRAMs for pedigree checks.
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Expected to generate the fingerprints file
        """
        return sample.get_cram_path().somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using a function from the `jobs` module.
        """
        cram_path = sample.get_cram_path()

        expected_path = self.expected_outputs(sample)
        j = somalier.extact_job(
            b=self.b,
            gvcf_or_cram_or_bam_path=cram_path,
            out_fpath=expected_path,
            overwrite=not self.check_intermediates,
            refs=self.refs,
            job_attrs=self.get_job_attrs(sample),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage(required_stages=CramSomalier)
class CramSomalierPedigree(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """

        prefix = dataset.get_bucket() / 'somalier' / dataset.alignment_inputs_hash()
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': dataset.get_web_bucket() / 'somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = {
            s.id: inputs.as_path(stage=CramSomalier, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.get_web_url():
            html_url = str(html_path).replace(str(dataset.get_web_bucket()), base_url)
        else:
            html_url = None

        jobs = somalier.pedigree(
            self.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.check_intermediates,
            out_samples_path=self.expected_outputs(dataset)['samples'],
            out_pairs_path=self.expected_outputs(dataset)['pairs'],
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=self.expected_outputs(dataset)['checks'],
            refs=self.refs,
        )
        return self.make_outputs(
            dataset, data=self.expected_outputs(dataset), jobs=jobs
        )


@stage(required_stages=CramSomalier)
class CramSomalierAncestry(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        *.somalier-ancestry.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """

        prefix = dataset.get_bucket() / 'ancestry' / dataset.alignment_inputs_hash()
        return {
            'tsv': prefix / f'{dataset.name}.somalier-ancestry.tsv',
            'html': dataset.get_web_bucket() / 'somalier-ancestry.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = {
            s.id: inputs.as_path(stage=CramSomalier, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.get_web_url():
            html_url = str(html_path).replace(str(dataset.get_web_bucket()), base_url)
        else:
            html_url = None

        j = somalier.ancestry(
            self.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.check_intermediates,
            out_tsv_path=self.expected_outputs(dataset)['tsv'],
            out_html_path=html_path,
            out_html_url=html_url,
            refs=self.refs,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[j])
