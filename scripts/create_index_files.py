#!/usr/bin/env python3

import click
from google.cloud import storage

from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_gcp_project, image_path
from cpg_utils.hail_batch import Batch, command, get_batch
from cpg_workflows.resources import STANDARD


def find_files_to_index(
    path: Path,
) -> list[tuple[str, str]]:
    """Finds bam and cram files and queues them for indexing if they are not already indexed."""
    client = storage.Client()
    bucket_name = path.as_uri().split('/')[2]
    prefix = '/'.join(path.as_uri().split('/')[3:])
    bucket = client.bucket(bucket_name, user_project=get_gcp_project())

    # Find all files with the given extensions
    files_to_index = []
    blobs = bucket.list_blobs(prefix=prefix)
    blob_names = [f'gs://{bucket_name}/{blob.name}' for blob in blobs]
    for blob_name in blob_names:
        # Check if the index file exists
        if blob_name.endswith('.bam'):
            index_blob_name = f'{blob_name}.bai'
        elif blob_name.endswith('.cram'):
            index_blob_name = f'{blob_name}.crai'
        else:
            continue
        if index_blob_name not in blob_names:
            files_to_index.append((blob_name, index_blob_name))
        else:
            continue

    return files_to_index


def index_with_samtools(
    b: Batch,
    input_files: list[tuple[str, str]],
    disk_size: str | None,
):
    """
    Index bam or cram files using samtools.
    """
    for file_to_index_path, index_file_path in input_files:
        input_file = b.read_input(file_to_index_path)
        j = b.new_bash_job(f'samtools index {file_to_index_path}')
        j.image(image_path('samtools'))
        # Set resource requirements
        if disk_size:
            storage_gb = int(disk_size)
        else:
            storage_gb = 20 if to_path(file_to_index_path).stat().st_size < 1e10 else 100
        nthreads = 8
        res = STANDARD.set_resources(
            j,
            ncpu=nthreads,
            storage_gb=storage_gb,
        )
        cmd = f'samtools index -@ {res.get_nthreads() - 1} {input_file} -o {j.output}'
        j.command(command(cmd, monitor_space=True))
        b.write_output(j.output, index_file_path)


@click.command()
@click.option(
    '--input-path',
    help='Path to the input files',
)
@click.option('--disk-size', help='Specify the disk size in gb')
def main(input_path: str, disk_size: str | None):
    # Find all bam and cram files that need to be indexed
    files_to_index = find_files_to_index(to_path(input_path))
    # Index the files
    index_with_samtools(get_batch(), files_to_index, disk_size)
    get_batch().run()


if __name__ == '__main__':
    main()
