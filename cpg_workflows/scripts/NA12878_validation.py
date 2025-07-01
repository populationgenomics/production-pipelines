"""
jobs relating to the validation steps of the pipeline
"""

from argparse import ArgumentParser

from cpg_utils import to_path, config
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import fasta_res_group, get_batch


def run_happy_on_gvcf(
    vcf_path: str,
    output_prefix: str,
):
    """Run hap.py on a single-sample gVCF (NA12878) using Truth data from config."""

    batch_instance = get_batch()

    # region: run BCFtools to filter the gVCF
    bcftools_job = batch_instance.new_job(f'Run BCFtools on {vcf_path!s}')
    bcftools_job.image(config.config_retrieve(['images', 'bcftools_120'])).memory('10Gi').storage('100Gi')

    # region: read input data into batch
    vcf_input = batch_instance.read_input_group(vcf=vcf_path, index=f'{vcf_path}.tbi')
    bcftools_job.declare_resource_group(
        vcf_output={
            'vcf': '{root}.vcf.gz',
            'index': '{root}.vcf.gz.tbi',
        },
    )
    # quick rinse of a gVCF to remove non-alt sites
    bcftools_job.command(f'bcftools view -c1 -Oz -o {bcftools_job.vcf_output["vcf"]} -W=tbi {vcf_input["vcf"]}')
    # endregion

    # region: run hap.py on the filtered VCF
    happy_j = batch_instance.new_job(f'Run Happy on {vcf_path!s}')
    happy_j.depends_on(bcftools_job)

    happy_j.image(config.config_retrieve(['images', 'hap-py'])).memory('100Gi').storage('100Gi').cpu(4)
    # read in sample-specific truth data from config
    truth_vcf = batch_instance.read_input(config_retrieve(['references', 'na12878','truth_vcf']))
    # declare all the expected outputs
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
    ref_genome_sdf = config.config_retrieve(['references', 'refgenome_sdf'])
    sdf = batch_instance.read_input_group(
        **{file.name: file.as_uri() for file in to_path(ref_genome_sdf).glob('*')},
    )

    # run the command
    happy_j.command(
        f"""
        mkdir {happy_j.output}
        
        hap.py {truth_vcf} {vcf_input["vcf"]}  \\
        -r {fasta_res_group(batch_instance)["base"]} \\ 
        -o {happy_j.output}/output \\
        --leftshift \\
        --threads 4 \\
        --preprocess-truth \ \
        --engine-vcfeval-path=/opt/hap.py/libexec/rtg-tools-install/rtg \\ 
        --engine-vcfeval-template {sdf} \\
        --engine=vcfeval 
        """
    )
    # endregion

    batch_instance.write_output(happy_j.output, output_prefix)

    return [bcftools_job, happy_j]


if __name__ == '__main__':
    parser = ArgumentParser(description='Run hap.py on NA12878 gVCF')
    parser.add_argument('vcf_path', type=str, help='Path to the gVCF file')
    parser.add_argument('output', type=str, help='Output prefix for the hap.py results')

    args = parser.parse_args()

    _jobs = run_happy_on_gvcf(vcf_path=args.vcf_path, output_prefix=args.output)

    get_batch().run(wait=False)
