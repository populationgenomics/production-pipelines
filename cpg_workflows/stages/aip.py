"""
This is a central script for the stages associated with AIP

AIP in general is supposed to be platform-agnostic, but this version will use
more explicitly CPG-coded alternatives for ease of execution

Takes as input:
    - Annotated MT path
    - Seqr - Internal mapping file
    - HPO.obo file
Generates:
    - PED file
    - Mapping between internal and external IDs
    - Latest participant panels
    - PanelApp results

Plan -
    - Generate PED file using Dataset method

Initial expectation - this will have to be run using Full permissions
as we will need to reference data in test

HPO file currently gs://cpg-common-test/references/aip/hpo_terms.obo
"""
from datetime import datetime
import logging
from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch


from cpg_workflows.workflow import (
    get_workflow,
    stage,
    Cohort,
    CohortStage,
    Dataset,
    DatasetStage,
    StageOutput,
    StageInput,
)
from cpg_workflows.jobs.aip.hpo_panel_match import main as panel_match_main


_GRANULAR_DATE: str | None = None
AIP_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_aip:latest'


def get_granular_date():
    """
    cached getter/setter for a simplified date representation
    """
    global _GRANULAR_DATE
    if _GRANULAR_DATE is None:
        _GRANULAR_DATE = get_config()['workflow'].get(
            'fake_date',
            datetime.now().strftime('%Y-%m-%d')
        )
    return _GRANULAR_DATE


@stage
class GeneratePanelData(DatasetStage):
    """
    PythonJob to find HPO-matched panels
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        only one output, the panel data
        """
        return {'hpo_panels': dataset.prefix() / 'reanalysis' / get_granular_date() / 'panel_data.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        py_job = get_batch().new_python_job('create_panel_data')
        py_job.image(AIP_IMAGE)

        expected_d = self.expected_outputs(dataset)
        hpo_file = get_batch().read_input('OBO FILE')  # todo
        py_job.call(
            panel_match_main,
            dataset.name,
            hpo_file,
            expected_d['hpo_panels'],
        )

        return self.make_outputs(dataset, data=expected_d, jobs=py_job)






