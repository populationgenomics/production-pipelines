"""
jobs relating to the validation steps of the pipeline
"""

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch


def get_sample_truth_data(sample_ext_id: str):
    """
    retrieve the reference truth specific to this individual

    Args:
        sample_ext_id (str):

    Returns:
        The specific truth data from config
    """

    ref_data = config_retrieve(['references', sample_ext_id])
    assert all(key in ref_data for key in ['vcf', 'index', 'bed'])
    return ref_data


def run_happy_on_vcf(
    vcf_path: str,
    sample_ext_id: str,
    out_prefix: str,
):
    """
    Run hap.py on the for this single-sample VCF
    Use Truth data from config

    Args:
        vcf_path (str): path to the single-sample VCF
        sample_ext_id (str): external ID to find reference data
        out_prefix (str): where to export happy outputs

    Returns:
        This Job or None
    """

    happy_j = get_batch().new_job(f'Run Happy on {sample_ext_id} VCF')
    happy_j.image(image_path('hap-py')).memory('100Gi').storage('100Gi').cpu(4)

    # region: read input data into batch
    vcf_input = get_batch().read_input_group(vcf=vcf_path, index=f'{vcf_path}.tbi')

    # read in sample-specific truth data from config
    ref_data = get_sample_truth_data(sample_ext_id)
    truth_input = get_batch().read_input_group(vcf=ref_data['vcf'], index=ref_data['index'], bed=ref_data['bed'])

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
    sdf = get_batch().read_input_group(
        **{file.name: file.as_uri() for file in to_path(ref_genome_sdf).glob('*')},
    )
    # endregion

    # includes the creation of the output temp file, which hap.py won't do
    command = (
        f'mkdir {happy_j.output} && '
        f'hap.py {truth_input["vcf"]} {vcf_input["vcf"]} '
        f'-r {fasta_res_group(get_batch())["base"]} -R {truth_input["bed"]} '
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
        batch_beds = get_batch().read_input_group(**strat_dict)
        command += f'--stratification {batch_beds["definition.tsv"]}'
    # endregion

    happy_j.command(command)

    get_batch().write_output(happy_j.output, out_prefix)
    return happy_j
