import click

from cpg_utils.hail_batch import Batch, command, get_batch, get_config, image_path


@click.command()
@click.option('--input_cram_path', required=True, type=str)
@click.option('--output_cram_path', required=True, type=str)
@click.option('--output_crai_path', required=True, type=str)
@click.option('--chr', required=True, type=str)
def main(input_cram_path: str, output_cram_path: str, output_crai_path: str, chr: str):
    """
    Subset a CRAM file to a single chromosome.
    """
    b = get_batch('subset_tob_cram')
    j = b.new_job('subset_tob_cram')
    j.image(image_path('samtools'))
    j.cpu(2)
    j.storage('150G')
    j.memory('standard')
    input_cram = b.read_input_group(**{'cram': input_cram_path, 'crai': input_cram_path + '.crai'})
    ref_path = b.read_input(
        'gs://cpg-common-main/references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta',
    )

    subset_cmd = f"""
    samtools view -T {ref_path} -C -o {j.output_cram} {input_cram['cram']} {chr} && \
    samtools index {j.output_cram} {j.output_cram}.crai
    """
    j.command(command(subset_cmd))
    b.write_output(j.output_cram, output_cram_path)
    b.write_output(j.output_cram + '.crai', output_crai_path)
    b.run(wait=False)


if __name__ == "__main__":
    main()
