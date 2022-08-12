#!/usr/bin/env python3

"""
Stages to run somalier tools.
"""

import logging

from cpg_utils.config import get_config

from .align import Align
from .. import Path
from ..targets import Dataset, Sample
from ..pipeline import stage, SampleStage, DatasetStage, StageInput, StageOutput
from ..jobs import somalier
from ..utils import exists

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
        return sample.make_cram_path().somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using a function from the `jobs` module.
        """
        cram_path = sample.make_cram_path()
        if get_config()['workflow'].get('check_inputs') and not cram_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(f'No CRAM found, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No CRAM found')

        expected_path = self.expected_outputs(sample)
        j = somalier.extact_job(
            b=self.b,
            gvcf_or_cram_or_bam_path=cram_path,
            out_somalier_path=expected_path,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage
class GvcfSomalier(SampleStage):
    """
    Genereate fingerprints from GVCFs for pedigree checks.
    """

    def expected_outputs(self, sample: Sample) -> Path:
        """
        Expected to generate the fingerprints file
        """
        return sample.make_gvcf_path().somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using a function from the `jobs` module.
        """
        gvcf_path = sample.make_gvcf_path()
        if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.warning(f'No GVCF found, skipping sample {sample}')
                return self.make_outputs(sample, skipped=True)
            else:
                return self.make_outputs(sample, error_msg=f'No GVCF found')

        expected_path = self.expected_outputs(sample)
        j = somalier.extact_job(
            b=self.b,
            gvcf_or_cram_or_bam_path=gvcf_path,
            out_somalier_path=expected_path,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            job_attrs=self.get_job_attrs(sample),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


# @stage
# class JointVcfSomalier(SampleStage):
#     """
#     Genereate fingerprints from joint VCF for pedigree checks.
#     """
#
#     def expected_outputs(self, cohort: Cohort) -> Path:
#         """
#         Expected to generate the fingerprints file
#         """
#         return sample.make_gvcf_path().somalier_path
#
#     def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
#         """
#         Using a function from the `jobs` module.
#         """
#         gvcf_path = sample.make_gvcf_path()
#         if get_config()['workflow'].get('check_inputs') and not gvcf_path.exists():
#             if get_config()['workflow'].get('skip_samples_with_missing_input'):
#                 logger.warning(f'No GVCF found, skipping sample {sample}')
#                 return self.make_outputs(sample, skipped=True)
#             else:
#                 return self.make_outputs(sample, error_msg=f'No GVCF found')
#
#         expected_path = self.expected_outputs(sample)
#         j = somalier.extact_job(
#             b=self.b,
#             gvcf_or_cram_or_bam_path=gvcf_path,
#             out_fpath=expected_path,
#             overwrite=not get_config()['workflow'].get('check_intermediates'),
#             job_attrs=self.get_job_attrs(sample),
#             sequencing_type=self.cohort.sequencing_type,
#         )
#         return self.make_outputs(sample, data=expected_path, jobs=[j])


EXCLUDE_HIGH_CONTAMINATION = False


@stage(required_stages=[Align])
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

        h = dataset.alignment_inputs_hash()
        prefix = dataset.prefix() / 'somalier' / 'cram' / h
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': dataset.web_prefix() / 'cram-somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Checks calls job from the pedigree module
        """
        verifybamid_by_sid = {}
        somalier_by_sid = {}
        for sample in dataset.get_samples():
            if EXCLUDE_HIGH_CONTAMINATION:
                verify_bamid_path = inputs.as_path(
                    stage=Align, target=sample, id='verify_bamid'
                )
                if not exists(verify_bamid_path):
                    logger.warning(
                        f'VerifyBAMID results {verify_bamid_path} do not exist for '
                        f'{sample}, somalier pedigree estimations might be affected'
                    )
                else:
                    verifybamid_by_sid[sample.id] = verify_bamid_path
            somalier_path = inputs.as_path(stage=Align, target=sample, id='somalier')
            if not exists(somalier_path):
                raise ValueError(
                    f'Somalier file does not exist for {sample}: {somalier_path}'
                )
            else:
                somalier_by_sid[sample.id] = somalier_path

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        if any(s.pedigree for s in dataset.get_samples()):
            expected_ped_path = dataset.write_ped_file(
                self.expected_outputs(dataset)['expected_ped']
            )
            jobs = somalier.pedigree(
                self.b,
                dataset,
                expected_ped_path=expected_ped_path,
                input_path_by_sid=somalier_by_sid,
                verifybamid_by_sid=verifybamid_by_sid,
                overwrite=not not get_config()['workflow'].get('check_intermediates'),
                out_samples_path=self.expected_outputs(dataset)['samples'],
                out_pairs_path=self.expected_outputs(dataset)['pairs'],
                out_html_path=html_path,
                out_html_url=html_url,
                out_checks_path=self.expected_outputs(dataset)['checks'],
                job_attrs=self.get_job_attrs(dataset),
                send_to_slack=True,
            )
            return self.make_outputs(
                dataset, data=self.expected_outputs(dataset), jobs=jobs
            )
        else:
            return self.make_outputs(dataset, skipped=True)


@stage(required_stages=[GvcfSomalier])
class GvcfSomalierPedigree(DatasetStage):
    """
    Checks pedigree from GVCF fingerprints.
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

        h = dataset.alignment_inputs_hash()
        prefix = dataset.prefix() / 'somalier' / 'gvcf' / h
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': dataset.web_prefix() / 'gvcf-somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = {
            s.id: inputs.as_path(stage=GvcfSomalier, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        if any(s.pedigree for s in dataset.get_samples()):
            expected_ped_path = dataset.write_ped_file(
                self.expected_outputs(dataset)['expected_ped']
            )
            jobs = somalier.pedigree(
                self.b,
                dataset,
                expected_ped_path=expected_ped_path,
                input_path_by_sid=fp_by_sid,
                overwrite=not not get_config()['workflow'].get('check_intermediates'),
                out_samples_path=self.expected_outputs(dataset)['samples'],
                out_pairs_path=self.expected_outputs(dataset)['pairs'],
                out_html_path=html_path,
                out_html_url=html_url,
                out_checks_path=self.expected_outputs(dataset)['checks'],
                job_attrs=self.get_job_attrs(dataset),
                send_to_slack=False,
            )
            return self.make_outputs(
                dataset, data=self.expected_outputs(dataset), jobs=jobs
            )
        else:
            return self.make_outputs(dataset, skipped=True)


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

        h = dataset.alignment_inputs_hash()
        prefix = dataset.prefix() / 'ancestry' / 'cram' / h
        return {
            'tsv': prefix / f'{dataset.name}.somalier-ancestry.tsv',
            'html': dataset.web_prefix() / 'somalier-ancestry.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = {
            s.id: inputs.as_path(stage=CramSomalier, target=s)
            for s in dataset.get_samples()
        }

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        j = somalier.ancestry(
            self.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            out_tsv_path=self.expected_outputs(dataset)['tsv'],
            out_html_path=html_path,
            out_html_url=html_url,
            job_attrs=self.get_job_attrs(dataset),
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[j])
