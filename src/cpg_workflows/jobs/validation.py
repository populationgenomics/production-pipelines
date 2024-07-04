"""
jobs relating to the validation steps of the pipeline
"""

import json
from csv import DictReader

from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_utils import to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import fasta_res_group, query_command
from cpg_workflows.workflow import SequencingGroup

from ..metamist import AnalysisStatus, get_metamist

SUMMARY_KEYS = {
    'TRUTH.TOTAL': 'true_variants',
    'METRIC.Recall': 'recall',
    'METRIC.Precision': 'precision',
}


def get_sample_truth_data(sequencing_group_id: str):
    """
    retrieve the reference truth specific to this individual

    Args:
        sequencing_group_id (str):

    Returns:
        The specific truth data from config
    """

    ref_data = get_config()['references'][sequencing_group_id]
    assert all(key in ref_data for key in ['vcf', 'index', 'bed'])
    return ref_data


def validation_mt_to_vcf_job(
    b: Batch,
    mt_path: str,
    sequencing_group_id: str,
    out_vcf_path: str,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None,
):
    """
    Take the single-dataset MT, and write to a VCF

    Args:
        b (hb.Batch): the batch to add jobs into
        mt_path (str): path of the AnnotateDataset MT
        sequencing_group_id (str): sequencing group name
        out_vcf_path (str): path to write new VCF to
        job_attrs (dict):
        depends_on (hb.Job|list[hb.Job]): jobs to depend on

    Returns:
        this single Job
    """
    from cpg_workflows.query_modules import validation

    vcf_j = b.new_job('VCF from dataset MT', (job_attrs or {}) | {'tool': 'hail query'})
    vcf_j.image(image_path('cpg_workflows'))
    vcf_j.command(
        query_command(
            validation,
            validation.single_sample_vcf_from_dataset_vcf.__name__,
            mt_path,
            sequencing_group_id,
            out_vcf_path,
            setup_gcp=True,
        ),
    )
    if depends_on:
        vcf_j.depends_on(*depends_on)

    return vcf_j


def run_happy_on_vcf(
    b: Batch,
    vcf_path: str,
    sequencing_group_ext_id: str,
    out_prefix: str,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None,
):
    """
    Run hap.py on the for this single-sample VCF
    Use Truth data from config

    Args:
        b (): the batch to create a new job in
        vcf_path (str): path to the single-sample VCF
        sequencing_group_ext_id (str): external ID to find reference data
        out_prefix (str): where to export happy outputs
        job_attrs ():
        depends_on ():

    Returns:
        This Job or None
    """

    happy_j = b.new_job(
        f'Run Happy on {sequencing_group_ext_id} VCF',
        (job_attrs or {}) | {'tool': 'hap.py'},
    )
    happy_j.image(image_path('hap-py')).memory('100Gi').storage('100Gi').cpu(4)
    if depends_on:
        happy_j.depends_on(*depends_on)

    # region: read input data into batch
    vcf_input = b.read_input_group(vcf=vcf_path, index=f'{vcf_path}.tbi')

    # read in sample-specific truth data from config
    ref_data = get_sample_truth_data(sequencing_group_ext_id)
    truth_input = b.read_input_group(vcf=ref_data['vcf'], index=ref_data['index'], bed=ref_data['bed'])

    happy_j.declare_resource_group(
        output={
            'happy_extended.csv': '{root}/output.extended.csv',
            'happy.vcf.gz': '{root}/output.vcf.gz',
            'happy.vcf.gz.tbi': '{root}/output.vcf.gz.tbi',
            'happy_roc.all.csv.gz': '{root}/output.roc.all.csv.gz',
            'happy_metrics.json.gz': '{root}/output.metrics.json.gz',
            'happy_runinfo.json': '{root}/output.runinfo.json',
            'summary.csv': '{root}/output.summary.csv',
        },
    )

    # read in all the SDF genome indices used by hap.py
    ref_genome_sdf = get_config()['references']['refgenome_sdf']
    sdf = b.read_input_group(
        **{file.name: file.as_uri() for file in to_path(ref_genome_sdf).glob('*')},
    )
    # endregion

    # includes the creation of the output temp file, which hap.py won't do
    command = (
        f'mkdir {happy_j.output} && '
        f'hap.py {truth_input["vcf"]} {vcf_input["vcf"]} '
        f'-r {fasta_res_group(b)["base"]} -R {truth_input["bed"]} '
        f'-o {happy_j.output}/output --leftshift '
        f'--threads 4 --preprocess-truth '
        f'--engine-vcfeval-path=/opt/hap.py/libexec/rtg-tools-install/rtg '
        f'--engine-vcfeval-template {sdf} --engine=vcfeval '
    )

    # region: stratification BED files
    if stratification := get_config()['references']['stratification']:
        strat_folder = to_path(stratification)
        strat_dict = {file.name: str(file) for file in strat_folder.glob('*')}
        assert 'definition.tsv' in strat_dict, 'definition.tsv file does not exist'
        batch_beds = b.read_input_group(**strat_dict)
        command += f'--stratification {batch_beds["definition.tsv"]}'
    # endregion

    happy_j.command(command)

    b.write_output(happy_j.output, out_prefix)
    return happy_j


def parse_and_post_results(
    vcf_path: str,
    sequencing_group_id: str,
    sequencing_group_ext_id: str,
    happy_csv: str,
    out_file: str,
):
    """
    Read the Hap.py results, and update Metamist
    This whole method is called as a scheduled PythonJob

    Args:
        vcf_path (str): path to the single-sample VCF
        sequencing_group_id (str): the SG ID
        sequencing_group_ext_id (str): the SG external ID
        happy_csv (str): CSV results from Hap.py
        out_file (str): where to write the JSON file

    Returns:
        this job
    """

    # handler for the CSV file
    happy_handle = to_path(happy_csv)

    ref_data = get_sample_truth_data(sequencing_group_id=sequencing_group_ext_id)

    # populate a dictionary of results for this sequencing group
    summary_data = {
        'type': 'validation_result',
        'query_vcf': vcf_path,
        'truth_vcf': ref_data['vcf'],
        'truth_bed': ref_data['bed'],
    }

    if stratification := get_config()['references']['stratification']:
        summary_data['stratified'] = stratification

    # read in the summary CSV file
    with happy_handle.open() as handle:
        summary_reader = DictReader(handle)
        for line in summary_reader:
            if line['Filter'] != 'PASS' or line['Subtype'] != '*':
                continue

            summary_key = f'{line["Type"]}_{line["Subset"]}'
            for sub_key, sub_value in SUMMARY_KEYS.items():
                summary_data[f'{summary_key}::{sub_value}'] = str(line[sub_key])

    with to_path(out_file).open('w', encoding='utf-8') as handle:
        json.dump(summary_data, handle)

    get_metamist().create_analysis(
        dataset=get_config()['workflow']['dataset'],
        status=AnalysisStatus('completed'),
        sequencing_group_ids=[sequencing_group_id],
        type_='qc',
        output=str(happy_handle.parent),
        meta=summary_data,
    )
