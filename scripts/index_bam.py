import click

import hail as hl

from cpg_utils.hail_batch import command, get_batch, image_path


@click.command()
@click.option('--bam-path', required=True, type=str)
@click.option('--output-path', required=True, type=str)
def main(bam_path: str, output_path: str):

    b = get_batch('Indexing bam')
    j = b.new_job('Index BAM')
    j.image(image_path('samtools'))
    j.storage('10Gi')
    bam_input = b.read_input(bam_path)
    cmd = f"samtools index {bam_input} {j.out_bam_index}"

    j.command(command(cmd))
    b.write_output(j.out_bam_index, output_path)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
