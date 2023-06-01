"""
jobs relating to the validation steps of the pipeline
"""


from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import fasta_res_group, image_path, query_command


def validation_mt_to_vcf_job(
    b: Batch,
    mt_path: Path,
    sample_id: str,
    out_vcf_path: Path,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None
):
    """
    Take the single-dataset MT, and write to a VCF

    Args:
        b (hb.Batch): the batch to add jobs into
        mt_path (str): path of the AnnotateDataset MT
        sample_id (str): sample name
        out_vcf_path (str): path to write new VCF to
        job_attrs (dict):
        depends_on (hb.Job|list[hb.Job]): jobs to depend on

    Returns:
        this single Job
    """
    from cpg_workflows.query_modules import validation

    vcf_j = b.new_job(
        f'VCF from dataset MT', (job_attrs or {}) | {'tool': 'hail query'}
    )
    vcf_j.image(image_path('cpg_workflows'))
    vcf_j.command(
        query_command(
            validation,
            validation.single_sample_vcf_from_dataset_vcf.__name__,
            str(mt_path),
            sample_id,
            out_vcf_path,
            setup_gcp=True,
        )
    )
    if depends_on:
        vcf_j.depends_on(*depends_on)

    return vcf_j


def run_happy_on_vcf(
    b: Batch,
    vcf_path: str,
    sample_ext_id: str,
    out_prefix: str,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None
):
    """
    Run hap.py on the for this single-sample VCF
    Use Truth data from config

    Args:
        b (): the batch to create a new job in
        vcf_path (str): path to the single-sample VCF
        sample_ext_id (str): external ID to find reference data
        out_prefix (str): where to export happy outputs
        job_attrs ():
        depends_on ():

    Returns:
        This Job or None
    """

    happy_j = b.new_job(
        f'Run Happy on {sample_ext_id} VCF', (job_attrs or {}) | {'tool': 'hap.py'}
    )
    happy_j.image(image_path('happy')).memory('100Gi').storage('100Gi').cpu(4)
    if depends_on:
        happy_j.depends_on(*depends_on)

    # region: read input data into batch
    vcf_input = b.read_input_group(vcf=vcf_path, index=f'{vcf_path}.tbi')

    ref_data = get_config()['references'][sample_ext_id]
    assert all(key in ref_data for key in ['vcf', 'index', 'bed'])

    # read in sample-specific truth data from config
    truth_input = b.read_input_group(
        vcf=ref_data['vcf'], index=ref_data['index'], bed=ref_data['bed']
    )

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
    if stratification := get_config()['stratification']:
        strat_folder = to_path(stratification)
        assert (
            strat_folder.exists()
        ), f'{stratification} does not exist, or was not accessible'

        definitions = strat_folder / 'definition.tsv'
        assert (
            definitions.exists()
        ), f'the region file {str(definitions)} does not exist'

        strat_bed_files = list(strat_folder.glob('*.bed*'))
        assert len(strat_bed_files) > 0, 'No bed files in the stratified BED folder'

        # create a dictionary to pass to input generation
        strat_dict = {'definition.tsv': str(definitions)}
        strat_dict.update({file.name: str(file) for file in strat_bed_files})
        batch_beds = b.read_input_group(**strat_dict)
        command += f'--stratification {batch_beds["definition.tsv"]}'
    # endregion

    happy_j.command(command)

    b.write_output(happy_j.output, out_prefix)
    return happy_j
