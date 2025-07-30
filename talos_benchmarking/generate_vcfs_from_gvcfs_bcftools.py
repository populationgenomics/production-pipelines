"""
Some gVCF files have been copied to test
Make them into VCFs
"""

import logging
import os

from cpg_utils import hail_batch, to_path
from cpg_workflows import utils


GS_GVCFS = 'gs://cpg-acute-care-test/talos_benchmarking/solo_gvcfs'
GS_VCFS = 'gs://cpg-acute-care-test/talos_benchmarking/solo_vcfs_bcftools'

# one batch to rule them all
batch_instance = hail_batch.get_batch()
# one ...reference to rule them all
reference = hail_batch.fasta_res_group(batch_instance)

bcftools_img = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.20-2'


def run_gvcf_to_vcf(input_sample: str, output_prefix: str):
    """
    Set up a job to run hap.py, generating a VCF from a GVCF.
    logic borrowed from https://github.com/populationgenomics/production-pipelines/blob/a93c4e5754930bff95f125b1292acbde95f7a777/cpg_workflows/scripts/NA12878_validation.py
    """

    sample_gvcf = f'{GS_GVCFS}/{input_sample}.g.vcf.gz'

    if not to_path(sample_gvcf).exists():
        logging.error(f'{sample_gvcf!s} does not exist, skipping')

    convert_job = batch_instance.new_job(f'Run pre.py on {sample_gvcf!s}')
    convert_job.image(bcftools_img).memory('10Gi').storage('10Gi')

    gvcf_input = batch_instance.read_input_group(gvcf=sample_gvcf, index=f'{sample_gvcf}.tbi')
    convert_job.declare_resource_group(
        vcf_output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        },
    )

    # use bcftools to convert the gVCF to a VCF
    # - normalise the gVCF to remove any multi-allelic sites (inc. splitting out variants from <NON_ALT> genotypes)
    # - annotate -x to drop the INFO fields (substantial in these gVCFs)
    # - convert to VCF using --gvcf2vcf
    # - compress the VCF and index it with tabix
    convert_job.command(
        f"""
        set -ex
        bcftools norm -m -any -f {reference["base"]} {gvcf_input["gvcf"]}  | \\
            bcftools annotate -x 'INFO' | \\
            bcftools view -c1 | \\
            bcftools convert --gvcf2vcf --fasta-ref {reference["base"]} \\
            -Oz -o {convert_job.vcf_output["vcf.gz"]} -W=tbi
        """
    )

    batch_instance.write_output(convert_job.vcf_output, output_prefix)


logging.basicConfig(level=logging.INFO)

# kick up a batch
hail_batch.init_batch()

# which gVCFs are available? Iterate over a generator
# use pre-py to convert gVCFs to VCFs
for each_gvcf in to_path(GS_GVCFS).glob('*.vcf.gz'):

    sample_id = each_gvcf.name.removesuffix('.g.vcf.gz')
    out_path = os.path.join(GS_VCFS, f'{sample_id}.vcf.gz')

    if utils.exists(out_path):
        logging.info(f'{out_path} exists, skipping')
        continue

    run_gvcf_to_vcf(input_sample=sample_id, output_prefix=os.path.join(GS_VCFS, sample_id))

hail_batch.get_batch().run(wait=False)
