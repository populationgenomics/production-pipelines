"""
This is a central script for the stages associated with Talos (previously AIP)

Due to privacy concerns, any configuration elements specific to an individual
project are stored in a separate config file, which is in a private repository.
 - production-pipelines-configuration/blob/main/config/talos/talos.toml

The cohort/project specific elements are described at the bottom of
the Talos default config file here:
 - production-pipelines/blob/main/configs/default/talos.toml

[cohorts.DATASET_NAME]
- cohort_panels, list[int], PanelApp panel IDs to apply to this project. By
  default, only the Mendeliome (137) is applied to all analyses.
- gene_prior, str, path to gene prior file to use for this project. When
  identifying Cat 2 (new disease-gene associations), the gene prior is used
  to determine whether a gene is novel or not. In ongoing analysis the gene
  prior used is the PanelApp result from the previous analysis.
- clinvar_filter, list[str], any clinvar submitters to filter out
  for this dataset (to blind this Talos run to its own clinvar entries).
  Note, this is not currently included in this workflow

[cohorts.DATASET_NAME.genome]  # repeated for exome if appropriate
- DATASET_NAME is taken from config[workflow][dataset]
- historic_results, str, path to historic results directory for this project. If
  present, new genes and new variant categories will be identified based on the
  state from the previous run. If absent, panel searching will default to using
  gene_prior above, and the results from this current run will not be stored in
  a suitable format to inform the next run.
- seqr_instance, str, URL for this seqr instance. Remove if not in seqr
- seqr_project, str, project ID for this project/seq type. Remove if not in seqr
- seqr_lookup, str, path to JSON file containing the mapping of CPG internal IDs
  to Seqr IDs, as generated by talos/helpers/process_seqr_metadata.py. If a new Seqr
  load has been performed, this file will need to be updated. This is stored in
  the test bucket so that it can be overwritten by a user. Remove if cohort is
  not in Seqr

Takes as input:
    - Annotated MT path (either found in metamist or directly in config)
    - HPO.obo file (read from cpg-common reference location)
    - Seqr<->Internal ID mapping file (if appropriate, from Config)
Generates:
    - PED file for this cohort (extended format, with additional columns for Ext. ID and HPO terms)
    - Latest panels found through HPO matches against this cohort
    - PanelApp results
    - Category-labelled VCF
    - Talos results JSON (metamist `aip-results` analysis)
    - Talos report HTML (metamist `aip-report` analysis)

This will have to be run using Full permissions as we will need to reference
data in test and main buckets.
"""

from datetime import datetime
from functools import lru_cache
from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import ConfigError, config_retrieve, image_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch
from cpg_workflows.resources import HIGHMEM, STANDARD
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import Dataset, DatasetStage, StageInput, StageOutput, stage
from metamist.graphql import gql, query

CHUNKY_DATE = datetime.now().strftime('%Y-%m-%d')
DATED_FOLDER = join('reanalysis', CHUNKY_DATE)
MTA_QUERY = gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)


@lru_cache(maxsize=None)
def query_for_sv_mt(dataset: str) -> list[tuple[str, str]]:
    """
    query for the latest SV MT for a dataset
    bonus - is we're searching for CNVs, we search for multiple
    return the full paths and filenames only, as 2 lists

    Args:
        dataset (str): project to query for

    Returns:
        str, the path to the latest MT for the given type
    """

    sv_type = 'cnv' if config_retrieve(['workflow', 'sequencing_type']) == 'exome' else 'sv'

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    # we only want the final stage MT, subset to the specific dataset
    final_stage_lookup = {'cnv': 'AnnotateDatasetCNV', 'sv': 'AnnotateDatasetSv'}

    result = query(MTA_QUERY, variables={'dataset': query_dataset, 'type': sv_type})
    mt_by_date: dict[str, str] = {}
    for analysis in result['project']['analyses']:
        if (
            analysis['output']
            and analysis['output'].endswith('.mt')
            and (analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type']))
            and (analysis['meta']['stage'] == final_stage_lookup[sv_type])
        ):
            mt_by_date[analysis['timestampCompleted']] = analysis['output']

    # perfectly acceptable to not have an input SV MT
    if not mt_by_date:
        return []

    if sv_type == 'cnv':
        full_paths = list(mt_by_date.values())
        filenames = [path.split('/')[-1] for path in full_paths]
        return list(zip(full_paths, filenames))

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    sv_file = mt_by_date[sorted(mt_by_date)[-1]]
    filename = sv_file.split('/')[-1]
    return [(sv_file, filename)]


@lru_cache(maxsize=None)
def query_for_latest_mt(dataset: str, entry_type: str = 'custom') -> str:
    """
    query for the latest MT for a dataset
    Args:
        dataset (str): project to query for
        entry_type (str): type of analysis entry to query for
    Returns:
        str, the path to the latest MT for the given type
    """

    # hot swapping to a string we can freely modify
    query_dataset = dataset

    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'
    result = query(MTA_QUERY, variables={'dataset': query_dataset, 'type': entry_type})
    mt_by_date = {}

    for analysis in result['project']['analyses']:
        if (
            analysis['output']
            and analysis['output'].endswith('.mt')
            and (analysis['meta']['sequencing_type'] == config_retrieve(['workflow', 'sequencing_type']))
        ):
            mt_by_date[analysis['timestampCompleted']] = analysis['output']

    if not mt_by_date:
        raise ValueError(f'No MT found for dataset {query_dataset}')

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return mt_by_date[sorted(mt_by_date)[-1]]


@lru_cache(2)
def get_clinvar_table(key: str = 'clinvar_decisions') -> str:
    """
    this is used to retrieve two types of object - clinvar_decisions & clinvar_pm5

    try and identify the clinvar table to use
    - try the config specified path
    - fall back to storage:common default path - try with multiple dir names in case we change this
    - if neither works, choose to fail instead

    Args
        key (str): the key to look for in the config

    Returns:
        a path to a clinvar table, or None
    """

    if (clinvar_table := config_retrieve(['workflow', key], None)) is not None:
        get_logger().info(f'Using clinvar table {clinvar_table} from config')
        return clinvar_table

    get_logger().info(f'No forced {key} table available, trying default')

    # try multiple path variations - legacy dir name is 'aip_clinvar', but this may also change
    for default_name in ['talos_clinvar', 'clinvarbitration', 'aip_clinvar']:
        clinvar_table = join(
            config_retrieve(['storage', 'common', 'analysis']),
            default_name,
            datetime.now().strftime('%y-%m'),
            f'{key}.ht',
        )

        if to_path(clinvar_table).exists():
            get_logger().info(f'Using clinvar table {clinvar_table}')
            return clinvar_table

    raise ValueError('no Clinvar Tables were identified')


@stage
class GeneratePED(DatasetStage):
    """
    this calls the script which reads pedigree data from metamist
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'pedigree': dataset.prefix() / DATED_FOLDER / 'pedigree.ped'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        generate a pedigree from metamist
        script to generate an extended pedigree format - additional columns for Ext. ID and HPO terms
        """
        job = get_batch().new_job('Generate PED from Metamist')
        job.cpu(0.25).memory('lowmem').image(image_path('talos'))
        expected_out = self.expected_outputs(dataset)
        query_dataset = dataset.name
        if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
            query_dataset += '-test'

        job.command(f'python3 talos/cpg_generate_pheno_ped.py {query_dataset} {job.output}')
        get_batch().write_output(job.output, str(expected_out["pedigree"]))
        get_logger().info(f'PED file for {dataset.name} written to {expected_out["pedigree"]}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=GeneratePED)
class GeneratePanelData(DatasetStage):
    """
    PythonJob to find HPO-matched panels
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        only one output, the panel data
        """
        return {'hpo_panels': dataset.prefix() / DATED_FOLDER / 'hpo_panel_data.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Find HPO-matched Panels: {dataset.name}')
        job.cpu(0.25).memory('lowmem').image(image_path('talos'))

        expected_out = self.expected_outputs(dataset)

        hpo_file = get_batch().read_input(config_retrieve(['workflow', 'obo_file']))
        local_ped = get_batch().read_input(str(inputs.as_path(target=dataset, stage=GeneratePED, key='pedigree')))
        job.command(
            'python3 talos/hpo_panel_match.py ' f'-i {local_ped!r} ' f'--hpo {hpo_file!r} ' f'--out {job.output}',
        )
        get_batch().write_output(job.output, str(expected_out["hpo_panels"]))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[GeneratePanelData])
class QueryPanelapp(DatasetStage):
    """
    query PanelApp for up-to-date gene lists
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'panel_data': dataset.prefix() / DATED_FOLDER / 'panelapp_data.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Query PanelApp: {dataset.name}')
        job.cpu(0.25).memory('lowmem').image(image_path('talos'))
        hpo_panel_json = inputs.as_path(target=dataset, stage=GeneratePanelData, key='hpo_panels')
        expected_out = self.expected_outputs(dataset)
        job.command(
            f'python3 talos/query_panelapp.py '
            f'--panels {get_batch().read_input(str(hpo_panel_json))!r} '
            f'--out_path {job.output} '
            f'--dataset {dataset.name} ',
        )
        get_batch().write_output(job.output, str(expected_out["panel_data"]))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED])
class RunHailFiltering(DatasetStage):
    """
    hail job to filter & label the MT

    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'labelled_vcf': dataset.prefix() / DATED_FOLDER / 'hail_labelled.vcf.bgz'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        input_mt = config_retrieve(['workflow', 'matrix_table'], query_for_latest_mt(dataset.name))
        job = get_batch().new_job(f'Run hail labelling: {dataset.name}')
        job.image(image_path('talos'))
        job.command('set -eux pipefail')

        # MTs can vary from <10GB for a small exome, to 170GB for a larger one
        # Genomes are more like 500GB
        seq_type: str = config_retrieve(['workflow', 'sequencing_type'], 'genome')
        required_storage: int = config_retrieve(['hail', 'storage', seq_type], 500)
        required_cpu: int = config_retrieve(['hail', 'cores', 'small_variants'], 8)
        HIGHMEM.set_resources(job, ncpu=required_cpu, storage_gb=required_storage)

        panelapp_json = get_batch().read_input(
            str(inputs.as_path(target=dataset, stage=QueryPanelapp, key='panel_data')),
        )
        pedigree = inputs.as_path(target=dataset, stage=GeneratePED, key='pedigree')
        expected_out = self.expected_outputs(dataset)

        # copy vcf & index out manually
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # peds can't read cloud paths
        local_ped = get_batch().read_input(str(pedigree))

        # find the clinvar tables, and localise
        clinvar_decisions = get_clinvar_table()
        clinvar_name = clinvar_decisions.split('/')[-1]

        # localise the clinvar decisions table
        job.command(
            f'cd $BATCH_TMPDIR && gcloud --no-user-output-enabled storage cp -r {clinvar_decisions} . && cd -',
        )

        # find, localise, and use the clinvar PM5 table
        pm5 = get_clinvar_table('clinvar_pm5')
        pm5_name = pm5.split('/')[-1]
        job.command(f'cd $BATCH_TMPDIR && gcloud --no-user-output-enabled storage cp -r {pm5} . && cd -')

        # finally, localise the whole MT (this takes the longest
        mt_name = input_mt.split('/')[-1]
        job.command(f'cd $BATCH_TMPDIR && gcloud --no-user-output-enabled storage cp -r {input_mt} . && cd -')

        job.command(
            f'python3 talos/hail_filter_and_label.py '
            f'--mt "${{BATCH_TMPDIR}}/{mt_name}" '
            f'--panelapp {panelapp_json!r} '
            f'--pedigree {local_ped!r} '
            f'--vcf_out {job.output["vcf.bgz"]} '
            f'--checkpoint "${{BATCH_TMPDIR}}/checkpoint.mt" '
            f'--clinvar "${{BATCH_TMPDIR}}/{clinvar_name}" '
            f'--pm5 "${{BATCH_TMPDIR}}/{pm5_name}" ',
        )
        get_batch().write_output(job.output, str(expected_out["labelled_vcf"]).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED])
class RunHailSVFiltering(DatasetStage):
    """
    hail job to filter & label the SV MT
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {
            filename: dataset.prefix() / DATED_FOLDER / f'label_{filename}.vcf.bgz'
            for _path, filename in query_for_sv_mt(dataset.name)
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        expected_out = self.expected_outputs(dataset)
        panelapp_json = get_batch().read_input(
            str(inputs.as_path(target=dataset, stage=QueryPanelapp, key='panel_data')),
        )
        pedigree = inputs.as_path(target=dataset, stage=GeneratePED, key='pedigree')

        # peddy can't read cloud paths
        local_ped = get_batch().read_input(str(pedigree))

        required_storage: int = config_retrieve(['hail', 'storage', 'sv'], 10)
        required_cpu: int = config_retrieve(['hail', 'cores', 'sv'], 2)
        sv_jobs: list = []
        for sv_path, sv_file in query_for_sv_mt(dataset.name):
            job = get_batch().new_job(f'Run hail SV labelling: {dataset.name}, {sv_file}')
            # manually extract the VCF and index
            job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
            job.image(image_path('talos'))
            STANDARD.set_resources(job, ncpu=required_cpu, storage_gb=required_storage, mem_gb=16)

            # copy the mt in
            job.command(f'gcloud --no-user-output-enabled storage cp -r {sv_path} .')
            job.command(
                f'python3 talos/hail_filter_sv.py '
                f'--mt {sv_file!r} '
                f'--panelapp {panelapp_json!r} '
                f'--pedigree {local_ped!r} '
                f'--vcf_out {job.output["vcf.bgz"]} ',
            )
            get_batch().write_output(job.output, str(expected_out[sv_file]).removesuffix('.vcf.bgz'))
            sv_jobs.append(job)

        return self.make_outputs(dataset, data=expected_out, jobs=sv_jobs)


@stage(
    required_stages=[GeneratePED, GeneratePanelData, QueryPanelapp, RunHailFiltering, RunHailSVFiltering],
    analysis_type='aip-results',  # note - legacy analysis type
    analysis_keys=['summary_json'],
)
class ValidateMOI(DatasetStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {'summary_json': dataset.prefix() / DATED_FOLDER / 'summary_output.json'}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Talos summary: {dataset.name}')
        job.cpu(2.0).memory('highmem').image(image_path('talos'))
        hpo_panels = get_batch().read_input(str(inputs.as_dict(dataset, GeneratePanelData)['hpo_panels']))
        pedigree = get_batch().read_input(str(inputs.as_path(target=dataset, stage=GeneratePED, key='pedigree')))
        hail_inputs = inputs.as_dict(dataset, RunHailFiltering)

        input_path = config_retrieve(['workflow', 'matrix_table'], query_for_latest_mt(dataset.name))

        # If there are SV VCFs, read each one in and add to the arguments
        sv_paths_or_empty = query_for_sv_mt(dataset.name)
        sv_vcf_arg = ''
        if sv_paths_or_empty:
            # only go looking for inputs from prior stage where we expect to find them
            hail_sv_inputs = inputs.as_dict(dataset, RunHailSVFiltering)
            for sv_path, sv_file in query_for_sv_mt(dataset.name):
                # bump input_path to contain both source files if appropriate
                # this is a string written into the metadata
                input_path += f', {sv_path}'

                labelled_sv_vcf = get_batch().read_input_group(
                    **{'vcf.bgz': str(hail_sv_inputs[sv_file]), 'vcf.bgz.tbi': f'{hail_sv_inputs[sv_file]}.tbi'},
                )['vcf.bgz']
                sv_vcf_arg += f' {labelled_sv_vcf} '

        if sv_vcf_arg:
            sv_vcf_arg = f'--labelled_sv {sv_vcf_arg}'

        panel_input = get_batch().read_input(str(inputs.as_dict(dataset, QueryPanelapp)['panel_data']))
        labelled_vcf = get_batch().read_input_group(
            **{'vcf.bgz': str(hail_inputs['labelled_vcf']), 'vcf.bgz.tbi': str(hail_inputs['labelled_vcf']) + '.tbi'},
        )['vcf.bgz']

        job.command(
            f'python3 talos/validate_categories.py '
            f'--labelled_vcf {labelled_vcf!r} '
            f'--out_json {job.output!r} '
            f'--panelapp {panel_input!r} '
            f'--pedigree {pedigree!r} '
            f'--input_path {input_path!r} '
            f'--participant_panels {hpo_panels!r} '
            f'--dataset {dataset.name!r} {sv_vcf_arg}',
        )
        get_batch().write_output(job.output, str(self.expected_outputs(dataset)['summary_json']))
        expected_out = self.expected_outputs(dataset)
        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[GeneratePED, ValidateMOI, QueryPanelapp, RunHailFiltering],
    analysis_type='aip-report',
    analysis_keys=['results_html', 'latest_html'],
    tolerate_missing_output=True,
)
class CreateTalosHTML(DatasetStage):
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {
            'results_html': dataset.prefix(category='web') / DATED_FOLDER / 'summary_output.html',
            'latest_html': dataset.prefix(category='web') / DATED_FOLDER / f'summary_latest_{CHUNKY_DATE}.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job(f'Talos HTML: {dataset.name}')
        job.cpu(1.0)
        job.memory('lowmem')
        job.image(image_path('talos'))

        results_json = get_batch().read_input(str(inputs.as_dict(dataset, ValidateMOI)['summary_json']))
        panel_input = get_batch().read_input(str(inputs.as_dict(dataset, QueryPanelapp)['panel_data']))
        expected_out = self.expected_outputs(dataset)

        # this will still try to write directly out - latest is optional, and splitting is arbitrary
        # Hail can't handle optional outputs being copied out AFAIK
        command_string = (
            f'python3 talos/html_builder.py '
            f'--results {results_json!r} '
            f'--panelapp {panel_input!r} '
            f'--output {str(expected_out["results_html"])!r} '
            f'--latest {str(expected_out["latest_html"])!r} '
            f'--dataset {dataset.name!r} '
        )

        if report_splitting := config_retrieve(['workflow', 'report_splitting', dataset.name], False):
            command_string += f' --split_samples {report_splitting}'

        job.command(command_string)

        # this + copy_common_env (called by default) should be enough to write using cloudpathlib
        authenticate_cloud_credentials_in_job(job)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


# probably shouldn't be recorded as a custom type
@stage(
    required_stages=ValidateMOI,
    analysis_keys=['seqr_file', 'seqr_pheno_file'],
    analysis_type='custom',
    tolerate_missing_output=True,
)
class GenerateSeqrFile(DatasetStage):
    """
    takes the results file from Seqr and produces a minimised form for Seqr ingestion
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        return {
            'seqr_file': dataset.prefix(category='analysis') / 'seqr_files' / f'{DATED_FOLDER}_seqr.json',
            'seqr_pheno_file': dataset.prefix(category='analysis') / 'seqr_files' / f'{DATED_FOLDER}_seqr_pheno.json',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        # pull out the config section relevant to this datatype & cohort

        # if it doesn't exist for this sequencing type, fail gracefully
        seq_type = config_retrieve(['workflow', 'sequencing_type'])
        try:
            seqr_lookup = config_retrieve(['cohorts', dataset.name, seq_type, 'seqr_lookup'])
        except ConfigError:
            get_logger().warning(f'No Seqr lookup file for {dataset.name} {seq_type}')
            return self.make_outputs(dataset, skipped=True)

        input_localised = get_batch().read_input(str(inputs.as_dict(dataset, ValidateMOI)['summary_json']))

        # create a job to run the minimisation script
        job = get_batch().new_job(f'Talos Prep for Seqr: {dataset.name}')
        job.cpu(1.0)
        job.memory('lowmem')
        job.image(image_path('talos'))
        lookup_in_batch = get_batch().read_input(seqr_lookup)
        job.command(
            f'python3 talos/minimise_output_for_seqr.py '
            f'{input_localised} {job.out_json} {job.pheno_json} '
            f'--external_map {lookup_in_batch}',
        )

        # write the results out
        expected_out = self.expected_outputs(dataset)
        get_batch().write_output(job.out_json, str(expected_out['seqr_file']))
        get_batch().write_output(job.pheno_json, str(expected_out['seqr_pheno_file']))
        return self.make_outputs(dataset, data=expected_out, jobs=job)