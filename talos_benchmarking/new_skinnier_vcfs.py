"""
pass each individual VCF through a my-variants-only filter
"""


from cpg_utils import hail_batch, to_path


batch_instance = hail_batch.get_batch('Slim down single sample VCFs')

GS_VCFS = 'gs://cpg-seqr-test/talos_benchmarking/solo_vcfs'

image = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/talos:PR_552'

output_dir = f'gs://cpg-acute-care-test/talos_benchmarking/new_slim_vcfs'

for input in to_path(GS_VCFS).glob('*.vcf.gz'):

    sample_id = input.name.split('.')[0]
    this_job = batch_instance.new_bash_job(f'Trim {sample_id}')
    this_job.image(image)
    this_job.cpu(1)

    vcf_input = batch_instance.read_input_group(gvcf=str(input), index=f'{input!s}.tbi').gvcf

    this_job.declare_resource_group(
        vcf_output={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        },
    )

    this_job.command(f'bcftools view -c1 -W=tbi -Oz -o {this_job.vcf_output["vcf.gz"]} {vcf_input}')

    batch_instance.write_output(this_job.vcf_output, f'{output_dir}/{sample_id}')

batch_instance.run(wait=False)
