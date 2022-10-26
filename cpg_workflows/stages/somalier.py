#!/usr/bin/env python3

"""
Stages to run somalier tools.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    DatasetStage,
    Dataset,
    exists,
)

from cpg_workflows.jobs import somalier
from .align import Align

EXCLUDE_HIGH_CONTAMINATION = False


@stage(required_stages=[Align])
class SomalierPedigree(DatasetStage):
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
        if get_config()['workflow'].get('skip_qc', False) is True:
            return {}
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
                    stage=Align, target=sample, key='verify_bamid'
                )
                if not exists(verify_bamid_path):
                    logging.warning(
                        f'VerifyBAMID results {verify_bamid_path} do not exist for '
                        f'{sample}, somalier pedigree estimations might be affected'
                    )
                else:
                    verifybamid_by_sid[sample.id] = verify_bamid_path
            somalier_path = inputs.as_path(stage=Align, target=sample, key='somalier')
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
