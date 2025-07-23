"""
Some gVCF files have been copied to test
Make them into VCFs
"""

import logging
import os

from cpg_utils import config, hail_batch, to_path
from cpg_workflows import utils


GS_GVCFS = 'gs://cpg-acute-care-test/talos_benchmarking/solo_gvcfs'
GS_VCFS = 'gs://cpg-acute-care-test/talos_benchmarking/solo_vcfs'

# one batch to rule them all
batch_instance = hail_batch.get_batch()
# one ...reference to rule them all
reference = hail_batch.fasta_res_group(batch_instance)


def run_gvcf_to_vcf(input_sample: str, output_prefix: str):
    """
    Set up a job to run hap.py, generating a VCF from a GVCF.
    logic borrowed from https://github.com/populationgenomics/production-pipelines/blob/a93c4e5754930bff95f125b1292acbde95f7a777/cpg_workflows/scripts/NA12878_validation.py
    """

    sample_gvcf = f'{GS_GVCFS}/{input_sample}.g.vcf.gz'

    if not to_path(sample_gvcf).exists():
        logging.error(f'{sample_gvcf!s} does not exist, skipping')

    prepy_j = batch_instance.new_job(f'Run pre.py on {sample_gvcf!s}')
    prepy_j.image(config.config_retrieve(['images', 'hap-py'])).memory('10Gi').storage('20Gi')


    gvcf_input = batch_instance.read_input_group(gvcf=sample_gvcf, index=f'{sample_gvcf}.tbi')
    prepy_j.declare_resource_group(
        vcf_output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        },
    )

    # convert the gVCF to a VCF, filtering to PASS sites and removing any variants genotyped as <NON_REF>
    prepy_j.command(
            f"""
            pre.py \\
            --convert-gvcf-to-vcf \\
            --filter-nonref \\
            --pass-only \\
            --reference {reference["base"]} \\
            {gvcf_input["gvcf"]} \\
            {prepy_j.vcf_output["vcf.gz"]}
            """,
    )

    batch_instance.write_output(prepy_j.vcf_output, output_prefix)


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