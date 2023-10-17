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
"""
from os.path import join
from datetime import datetime
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    image_path,
)
from metamist.graphql import gql, query

from cpg_workflows.batch import get_batch
from cpg_workflows.workflow import (
    stage,
    Dataset,
    DatasetStage,
    StageOutput,
    StageInput,
)
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.jobs.aip.hpo_panel_match import main as panel_match_main

RUN_CONFIG = get_config()
CHUNKY_DATE = datetime.now().strftime('%Y-%m-%d')
DATED_FOLDER = join('reanalysis', CHUNKY_DATE)
MTA_QUERY = gql(
    """
    query MyQuery($dataset: String!) {
        project(name: $dataset) {
            analyses(type: {eq: "custom"}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
            }
        }
    }
""")
# validate(MTA_QUERY)


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
        py_job.image(get_config()['workflow']['driver_image'])

        # only needs a really tiny image
        py_job.cpu(0.25).memory('lowmem')
        expected_d = self.expected_outputs(dataset)
        hpo_file = get_batch().read_input(get_config()['workflow']['obo_file'])
        py_job.call(
            panel_match_main, dataset.name, hpo_file, str(expected_d['hpo_panels'])
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
        job.cpu(0.25).memory('lowmem')
        job.image(image_path('aip'))

        # auth and copy env
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)

        hpo_panel_json = inputs.as_path(
            target=dataset, stage=GeneratePanelData, key='hpo_panels'
        )
        expected_out = self.expected_outputs(dataset)
        job.command(
            f'python3 reanalysis/query_panelapp.py '
            f'--panels "{str(hpo_panel_json)}" '
            f'--out_path "{str(expected_out["panel_data"])}" '
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
    result = query(MTA_QUERY, variables={'dataset': dataset})
    mt_by_date = {}
    seq_type_exome = get_config()['workflow'].get('sequencing_type') == 'exome'
    for analysis in result['project']['analyses']:
        if analysis['output'] and analysis['output'].endswith('.mt') and (
            (seq_type_exome and '/exome/' in analysis['output'])
            or (not seq_type_exome and '/exome/' not in analysis['output'])
        ):
            mt_by_date[analysis['timestampCompleted']] = analysis['output']

    if not mt_by_date:
        raise ValueError(f'No MT found for dataset {dataset}')

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return mt_by_date[sorted(mt_by_date)[-1]]


@stage(required_stages=[QueryPanelapp])
class RunHailFiltering(DatasetStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {
            'labelled_vcf': dataset.prefix() / DATED_FOLDER / 'hail_labelled.vcf.bgz',
            'pedigree': dataset.prefix() / DATED_FOLDER / 'pedigree.ped',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        # either do as you're told, or find the latest
        # this could potentially follow on from the AnnotateDataset stage
        # if this becomes integrated in the main pipeline
        input_mt = get_config()['workflow'].get(
            'matrix_table', query_for_latest_mt(dataset.name)
        )
        job = get_batch().new_job('run hail labelling')
        job.image(image_path('aip'))
        HIGHMEM.set_resources(job, mem_gb=16.0)

        # auth and copy env
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)

        panelapp_json = inputs.as_path(
            target=dataset, stage=QueryPanelapp, key='panel_data'
        )
        expected_out = self.expected_outputs(dataset)
        pedigree = dataset.write_ped_file(out_path=expected_out['pedigree'])
        # peddy can't read cloud paths
        local_ped = get_batch().read_input(str(pedigree))
        job.command(
            f'python3 reanalysis/hail_filter_and_label.py '
            f'--mt "{input_mt}" '
            f'--panelapp "{panelapp_json}" '
            f'--pedigree "{local_ped}" '
            f'--vcf_out "{str(expected_out["labelled_vcf"])}" '
        )

        return self.make_outputs(dataset, data=expected_out, jobs=job)


def _aip_summary_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Add meta.type to custom analysis object
    """
    return {'type': 'aip_output_json'}


@stage(
    required_stages=[GeneratePanelData, QueryPanelapp, RunHailFiltering],
    analysis_type='custom',
    analysis_keys=['summary_json'],
    update_analysis_meta=_aip_summary_meta,
)
class ValidateMOI(DatasetStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'summary_json': dataset.prefix() / DATED_FOLDER / 'summary_output.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'AIP summary for {dataset.name}')
        job.cpu(1.0).memory('standard')

        # auth and copy env
        job.image(image_path('aip'))
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)
        hpo_panels = str(inputs.as_dict(dataset, GeneratePanelData)['hpo_panels'])
        hail_inputs = inputs.as_dict(dataset, RunHailFiltering)
        panel_input = str(inputs.as_dict(dataset, QueryPanelapp)['panel_data'])
        # peddy can't read cloud paths
        local_ped = get_batch().read_input(str(hail_inputs['pedigree']))
        labelled_vcf = str(hail_inputs['labelled_vcf'])
        out_json_path = str(self.expected_outputs(dataset)['summary_json'])
        input_mt = get_config()['workflow'].get(
            'matrix_table', query_for_latest_mt(dataset.name)
        )
        job.command(
            f'python3 reanalysis/validate_categories.py '
            f'--labelled_vcf "{labelled_vcf}" '
            f'--out_json "{out_json_path}" '
            f'--panelapp "{panel_input}" '
            f'--pedigree "{local_ped}" '
            f'--input_path "{input_mt}" '
            f'--participant_panels "{hpo_panels}" '
        )
        expected_out = self.expected_outputs(dataset)
        return self.make_outputs(dataset, data=expected_out, jobs=job)


def _aip_html_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, str]:
    """
    Add meta.type to custom analysis object
    This isn't quite conformant with what AIP alone produces
    e.g. doesn't have the full URL to the results in GCP
    """
    return {
        'type': 'aip_output_json',
        'is_singleton': '',  # not doing this at the moment
        'is_exome': get_config()['workflow'].get('sequencing_type') == 'exome'
    }


@stage(
    required_stages=[ValidateMOI],
    analysis_type='custom',
    analysis_keys=['results_html', 'latest_html'],
    update_analysis_meta=_aip_summary_meta,
)
class CreateAIPHTML(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {
            'results_html': dataset.prefix(category='web') / DATED_FOLDER / 'summary_output.html',
            'latest_html': dataset.prefix(category='web') / DATED_FOLDER / f'summary_latest_{CHUNKY_DATE}.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'AIP HTML for {dataset.name}')
        job.cpu(1.0).memory('lowmem')

        # auth and copy env
        authenticate_cloud_credentials_in_job(job)
        copy_common_env(job)

        moi_inputs = inputs.as_dict(dataset, ValidateMOI)['summary_json']
        panel_input = inputs.as_dict(dataset, QueryPanelapp)['panel_data']
        pedigree = inputs.as_dict(dataset, RunHailFiltering)['pedigree']
        local_ped = get_batch().read_input(str(pedigree))

        expected_out = self.expected_outputs(dataset)
        job.command(
            f'python3 reanalysis/html_builder.py '
            f'--results "{moi_inputs}" '
            f'--panelapp "{panel_input}" '
            f'--pedigree "{local_ped}" '
            f'--output "{expected_out["results_html"]}" '
            f'--latest "{expected_out["latest_html"]}" '
        )

        expected_out = self.expected_outputs(dataset)
        return self.make_outputs(dataset, data=expected_out, jobs=job)
