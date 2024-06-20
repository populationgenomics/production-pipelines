import click

import hail as hl

from cpg_utils.hail_batch import command, get_batch, image_path


@click.command()
@click.option('--gvcf-path', required=True, type=str)
@click.option('--output-path', required=True, type=str)
def main(gvcf_path: str, output_path: str):

    b = get_batch('Subsample VCF')
    j = b.new_job('Subsample VCF')
    j.image(image_path('bcftools'))
    j.storage('10Gi')
    gvcf_input = b.read_input_group(gvcf=gvcf_path, gvcf_index=f'{gvcf_path}.tbi')

    cmd = f"bcftools view -r chr2:25227855-25342590 {gvcf_input.gvcf} -Oz -o {j.out_vcf}"

    j.command(command(cmd))
    b.write_output(j.out_vcf, output_path)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
