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
from os.path import join
from datetime import datetime
from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, copy_common_env
from metamist.graphql import gql, query
from cpg_workflows.batch import get_batch

from cpg_workflows.workflow import (
    stage,
    Dataset,
    DatasetStage,
    StageOutput,
    StageInput,
)
from cpg_workflows.jobs.aip.hpo_panel_match import main as panel_match_main


DATED_FOLDER = join('reanalysis', get_config()['workflow'].get(
    'fake_date',
    datetime.now().strftime('%Y-%m-%d')
))

# TODO move to default config
AIP_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_aip:latest'


@stage
class GeneratePanelData(DatasetStage):
    """
    PythonJob to find HPO-matched panels
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        only one output, the panel data
        """
        return {'hpo_panels': dataset.prefix() / DATED_FOLDER / 'hpo_panel_data.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        py_job = get_batch().new_python_job('create_panel_data')
        py_job.image(AIP_IMAGE)

        expected_d = self.expected_outputs(dataset)
        hpo_file = get_batch().read_input('OBO FILE')  # todo
        py_job.call(
            panel_match_main,
            dataset.name,
            hpo_file,
            expected_d['hpo_panels']
        )

        return self.make_outputs(dataset, data=expected_d, jobs=py_job)


@stage(required_stages=[GeneratePanelData])
class QueryPanelapp(DatasetStage):
    """
    query PanelApp for up-to-date gene lists
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'panel_data': dataset.prefix() / DATED_FOLDER / 'panelapp_data.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        job = get_batch().new_job('query panel data')
        job.image(AIP_IMAGE)

        # auth and copy env
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)

        hpo_panel_json = inputs.as_str(target=dataset, stage=GeneratePanelData, key='hpo_panels')
        expected_out = self.expected_outputs(dataset)
        job.command(
            f'python3 reanalysis/query_panelapp.py '
            f'--panels {hpo_panel_json} '
            f'--out_json {expected_out["panel_data"]}'
        )

        return self.make_outputs(dataset, data=expected_out, jobs=job)


def query_for_latest_mt(dataset: str) -> str:
    """
    query for the latest MT for a dataset

    Args:
        dataset (str):

    Returns:
        str, the path to the latest MT
    """
    my_query = gql(
        """
        query MyQuery($dataset: String!) {
            project(name: $dataset) {
                analyses(type: {eq: "custom"}, status: {eq: COMPLETED}) {
                    output
                }
            }
        }
    """)
    result = query(my_query, variables={'dataset': dataset})
    seq_type_exome = get_config()['workflow'].get('sequencing_type') == 'exome'
    mt_path = ''
    for analysis in result['project']['analyses']:
        if (analysis['output'].endswith('.mt') and
                (
                        (seq_type_exome and '/exome/' in analysis['output']) or
                        (not seq_type_exome and '/exome/' not in analysis['output'])
                )
        ):
            mt_path = analysis['output']

    if not mt_path:
        raise ValueError(f'No MT found for dataset {dataset}')
    return mt_path


@stage(required_stages=[QueryPanelapp])
class RunHailFiltering(DatasetStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'labelled_vcf': dataset.prefix() / DATED_FOLDER / 'hail_labelled.vcf.bgz',
                'pedigree': dataset.prefix() / DATED_FOLDER / 'pedigree.ped'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        # either do as you're told, or find the latest
        # this could potentially follow on from the AnnotateDataset stage
        input_mt = get_config()['workflow'].get('matrix_table', query_for_latest_mt(dataset.name))
        job = get_batch().new_job('run hail labelling')
        job.image(AIP_IMAGE)
        job.memory('32Gi')

        # auth and copy env
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)

        panelapp_json = inputs.as_str(target=dataset, stage=QueryPanelapp, key='panel_data')
        expected_out = self.expected_outputs(dataset)
        pedigree = dataset.write_ped_file(out_path=expected_out['pedigree'])
        local_ped = get_batch().read_input(pedigree)  # peddy can't read cloud paths
        job.command(
            f'python3 reanalysis/hail_filter_and_label.py '
            f'--mt {input_mt} '
            f'--panelapp {panelapp_json} '
            f'--pedigree {local_ped} '
            f'--vcf_out {expected_out["labelled_vcf"]}'
        )

        return self.make_outputs(dataset, data=expected_out, jobs=job)
