"""
takes the ms VCFs we have available, and trims them
initial benchmarking run showed it was faster to run on single sample VCFs, including merging
that's probably contaminated by variants where none of the samples have a positive call, and these
are artefacts from the MatrixTable to VCF export

Evidence for this:

- I didn't intentionally stop it happening
- the 5-sample MT contains as many variants as the 250-sample MT
- the MT VCF in the single-sample merge pipeline contains 60% of the samples in the pre-merged VCF pipeline

So...

read through the multisample VCFs, and trim them to remove INFO fields (for overall data size reduction)
and remove any variants where the contained samples have no positive calls
"""

import logging

from cpg_utils import hail_batch, to_path


logging.basicConfig(level=logging.INFO)

bcftools_img = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.20-2'
for each_group in [5, 10, 25, 50, 100, 250]:

    output_folder = f'gs://cpg-acute-care-test/talos_benchmarking/ms_vcfs_trimmed/{each_group}'

    if to_path(f'{output_folder}/{each_group}.vcf.gz').exists():
        logging.info(f'{output_folder}/{each_group}.vcf.gz exists, skipping')
        continue

    ms_vcf = f'gs://cpg-acute-care-test/talos_benchmarking/ms_vcfs/{each_group}.vcf.bgz'
    local_vcf = hail_batch.get_batch().read_input(ms_vcf)

    new_job = hail_batch.get_batch().new_bash_job(f'Run BCFtools trimming for {each_group} MS VCF')
    new_job.cpu(16).memory('32GiB').storage('250GiB')
    new_job.image(bcftools_img)

    new_job.declare_resource_group(
        vcf_output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        },
    )

    new_job.command(
        f"""
    set -x

    bcftools annotate -x 'INFO' {local_vcf} | \\
        bcftools view -c1 -Oz -o {new_job.vcf_output["vcf.gz"]} -W=tbi
    """
    )

    hail_batch.get_batch().write_output(new_job.vcf_output, f'{output_folder}/{each_group}')

hail_batch.get_batch().run(wait=False)
