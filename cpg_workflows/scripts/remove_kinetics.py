import click
from google.cloud import storage

from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch, init_batch


@click.command()
@click.option('--bucket', type=str, required=True, help='GCS bucket name where input ubams are stored.')
@click.option('--path-prefix', type=str, required=True, help='Path prefix within the bucket to look for ubams.')
@click.option('--output-dir', required=True, help='Path to save the output dataset')
def main(bucket: str, path_prefix: str, output_dir: str):

    # Get list of ubams from input_dir in GCS
    client = storage.Client()
    bucket_name = bucket
    prefix = path_prefix

    blobs = client.list_blobs(bucket_name, prefix=prefix)
    ubam_files: list[str] = []
    for blob in blobs:
        name = blob.name
        # ubams have the same file extension as a regular bam
        # We deliberately skip files with no_kinetics.bam suffix to avoid re-processing
        if name.endswith('no_kinetics.bam'):
            continue
        if name.endswith('.bam'):
            ubam_files.append(f'gs://{bucket_name}/{name}')
    print(f"Found {ubam_files} bams")

    # Initialise Hail Batch
    init_batch(
        worker_memory=config_retrieve(['workflow', 'worker_memory'], 'highmem'),
        driver_memory=config_retrieve(['workflow', 'driver_memory'], 'highmem'),
        driver_cores=config_retrieve(['workflow', 'driver_cores'], 2),
    )

    b = get_batch()

    # Create job for each ubam to remove kinetics data
    for ubam_path in ubam_files:
        j = b.new_bash_job(f'Remove kinetics data for {ubam_path}')
        j.image(image_path('samtools'))
        j.cpu(4)
        j.storage('600Gi')
        j.memory('16Gi')

        input_bam = b.read_input(ubam_path)
        j.command(
            f"""
            samtools view -@ 4 --remove-tag=fi,fn,fp,ri,rn,rp -b {input_bam} -o {j.output}
            """,
        )

        print(f'Writing output to {output_dir}/{ubam_path.split("/")[-1]}'.replace('.bam', '.no_kinetics.bam'))
        b.write_output(j.output, f'{output_dir}/{ubam_path.split("/")[-1]}'.replace('.bam', '.no_kinetics.bam'))

    b.run(wait=True)


if __name__ == '__main__':
    main()
