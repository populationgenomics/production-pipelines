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
) -> list[Path]:
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

        index_blob = bucket.blob(index_blob_name)
        if not index_blob.exists():
            files_to_index.append(to_path(blob_name))
        else:
            print(f'Index file already exists: {index_blob_name}')

    return files_to_index


def index_files(
    b: Batch,
    files_to_index: list[Path],
):
    """
    Index a list of bam or cram files using samtools.
    """
    for file_to_index in files_to_index:
        if file_to_index.suffix == '.bam':
            f = b.read_input_group(bam=file_to_index)
            index_with_samtools(b, f)
        elif file_to_index.suffix == '.cram':
            f = b.read_input_group(cram=file_to_index)
            index_with_samtools(b, f)
        else:
            raise ValueError(f'Unknown file extension: {file_to_index.suffix}')


def index_with_samtools(
    b: Batch,
    file_to_index: ResourceGroup,
):
    """
    Index a bam or cram file using samtools and return the job and the index file.
    """
    assert isinstance(file_to_index, ResourceGroup)

    job_name = 'Samtools index'
    j_attrs = dict(label=job_name, tool='samtools')
    j = b.new_job(name=job_name, attributes=j_attrs)
    j.image(image_path('samtools'))

    # Set resource requirements
    nthreads = 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=100,
    )
    if file_to_index.bam:
        j.declare_resource_group(
            out_bam={
                'bam': '{root}.bam',
                'bam.bai': '{root}.bam.bai',
            },
        )
        cmd = f'samtools index -@ {res.get_nthreads() - 1} {j.bam} -o {j.out_bam["bam.bai"]}'
        j.command(command(cmd, monitor_space=True))
        b.write_output(j.out_bam, f'{j.out_bam["bam.bai"].removesuffix(".bai")}')
    elif file_to_index.cram:
        j.declare_resource_group(
            out_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            },
        )
        cmd = f'samtools index -@ {res.get_nthreads() - 1} {j.cram} -o {j.out_cram["cram.crai"]}'
        j.command(command(cmd, monitor_space=True))
        b.write_output(j.out_cram, f'{j.out_cram["cram.crai"].removesuffix(".crai")}')
    else:
        raise ValueError('Resource group must contain a bam or cram file')


@click.command()
@click.option(
    '--input-path',
    help='Path to the input files',
)
def main(input_path: str):
    # Find all bam and cram files that need to be indexed
    files_to_index = find_files_to_index(to_path(input_path))
    # Index the files
    index_files(get_batch(), files_to_index)
    get_batch().run()


if __name__ == '__main__':
    main()
