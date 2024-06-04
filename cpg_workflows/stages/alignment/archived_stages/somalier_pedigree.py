import logging

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import somalier
from cpg_workflows.stages.alignment.cram_qc import CramQC
from cpg_workflows.targets import Dataset
from cpg_workflows.utils import exists
from cpg_workflows.workflow import (
    DatasetStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(required_stages=[CramQC])
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
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        prefix = dataset.prefix() / 'somalier' / 'cram' / dataset.alignment_inputs_hash()
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
        verifybamid_by_sgid = {}
        somalier_path_by_sgid = {}
        for sequencing_group in dataset.get_sequencing_groups():
            if config_retrieve(['somalier', 'exclude_high_contamination'], False):
                verify_bamid_path = inputs.as_path(stage=CramQC, target=sequencing_group, key='verify_bamid')
                if not exists(verify_bamid_path):
                    logging.warning(
                        f'VerifyBAMID results {verify_bamid_path} do not exist for '
                        f'{sequencing_group}, somalier pedigree estimations might be affected',
                    )
                else:
                    verifybamid_by_sgid[sequencing_group.id] = verify_bamid_path
            somalier_path = inputs.as_path(stage=CramQC, target=sequencing_group, key='somalier')
            somalier_path_by_sgid[sequencing_group.id] = somalier_path

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        if any(sg.pedigree.dad or sg.pedigree.mom for sg in dataset.get_sequencing_groups()):
            expected_ped_path = dataset.write_ped_file(self.expected_outputs(dataset)['expected_ped'])
            jobs = somalier.pedigree(
                b=get_batch(),
                dataset=dataset,
                expected_ped_path=expected_ped_path,
                somalier_path_by_sgid=somalier_path_by_sgid,
                verifybamid_by_sgid=verifybamid_by_sgid,
                out_samples_path=self.expected_outputs(dataset)['samples'],
                out_pairs_path=self.expected_outputs(dataset)['pairs'],
                out_html_path=html_path,
                out_html_url=html_url,
                out_checks_path=self.expected_outputs(dataset)['checks'],
                job_attrs=self.get_job_attrs(dataset),
                send_to_slack=True,
            )
            return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
        else:
            return self.make_outputs(dataset, skipped=True)
