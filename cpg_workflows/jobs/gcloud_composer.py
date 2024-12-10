"""
All jobs required to take a shard manifest (containing paths to VCF fragments) and produce a single VCF file

Gcloud Storage Compose is a can concatenate GCP objects without localising, and is limited to 32 objects per command
There is a limitation here that Condense can't run across buckets. Here we address this by keeping fragments and
intermediates in a temp bucket, then moving the final output to the desired location.

This implementation is designed to take a manifest of VCF fragments, and produce a single VCF file
Though conceptually this can be generalised to 'take a manifest of objects, and produce a single object'

The implementation is as follows:
1. Read the manifest file
2. Split the manifest into chunks of 32 objects
3. For each chunk, run a gcloud compose command to concatenate the objects into an intermediate
4. Keep condensing the intermediates until there is only one object left
5. Run a final job bcftools command to view, compress, and index the resulting VCF
6. Write the final VCF to a permanent GCS bucket
"""

from os.path import join

import hailtop.batch as hb

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, image_path
from cpg_workflows.utils import chunks, get_logger


def get_gcloud_condense_command_string(fragment_list: list[str], output: str) -> str:
    """
    Generate a gcloud compose command string to concatenate a list of VCF fragments
    write the product to the output path
    Args:
        fragment_list ():
        output ():
    Returns:
    """

    cat_string = 'gcloud storage objects compose '
    for fragment in fragment_list:
        cat_string += f'{fragment} '
    cat_string += f'{output}'
    return cat_string


def compose_condense_fragments(
    path_list: list[str],
    temp_prefix: str,
    chunk_size: int = config_retrieve(['gcloud_condense', 'chunk_size'], 32),
    job_attrs: dict | None = None,
) -> tuple[list[str], hb.batch.job.Job]:
    """
    takes a list of things to condense, creates jobs to condense them. Returns a list of the result paths

    Args:
        path_list (list[str]): all the paths to condense
        temp_prefix ():
        chunk_size (): number of objects to condense at once
        job_attrs (dict): job attributes
    Returns:
        list of paths to the condensed objects
    """

    job_name = job_attrs.get('name', 'ChunkBuster') if job_attrs else 'ChunkBuster'
    chunk_job = get_batch().new_job(name=f'{job_name}_{len(path_list)}', attributes=job_attrs)
    chunk_job.image(config_retrieve(['workflow', 'driver_image']))
    authenticate_cloud_credentials_in_job(chunk_job)

    sub_chunks = []
    for i, chunk in enumerate(chunks(path_list, chunk_size=chunk_size)):
        # either the named file, or a temp location
        chunk_output = join(temp_prefix, f'chunk_{i}.vcf.bgz')
        chunk_job.command(get_gcloud_condense_command_string(chunk, chunk_output))
        sub_chunks.append(chunk_output)
    return sub_chunks, chunk_job


def gcloud_compose_vcf_from_manifest(
    manifest_path: Path,
    intermediates_path: str,
    output_path: str,
    job_attrs: dict,
) -> list[hb.batch.job.Job]:
    """
    compose a series gcloud commands to condense a list of VCF fragments into a single VCF
    gcloud compose has a limit on 32 separate inputs, so at higher input counts we need to do a rolling merge

    compose can only run within a bucket, so we can't do all the intermediate generation in one bucket, and
    then move to a permanent location

    Args:
        manifest_path (Path): path to a file in GCP, one GCP file path per line
        intermediates_path (str): path to a bucket, where intermediate files will be written
        output_path (str): path to write the final file to. Must be in the same bucket as the manifest
        job_attrs (dict): job attributes

    Returns:
        a list of the jobs which will generate a single output from the manifest of fragment paths
    """

    get_logger().info(f'Reading manifest from {manifest_path}')

    # prefix these shard names to get full GCP path for each
    fragment_files = [
        join(str(manifest_path.parent), fragment.strip()) for fragment in manifest_path.open().readlines()
    ]

    # rolling squash of the chunks, should enable infinite-ish scaling
    temp_chunk_prefix_num = 1
    condense_jobs: list[hb.batch.job.Job] = []
    while len(fragment_files) > 1:
        condense_temp = join(intermediates_path, f'temp_chunk_{temp_chunk_prefix_num}')
        temp_chunk_prefix_num += 1
        fragment_files, condense_job = compose_condense_fragments(
            path_list=fragment_files,
            temp_prefix=condense_temp,
            job_attrs=job_attrs,
        )
        if condense_jobs:
            condense_job.depends_on(condense_jobs[-1])

        condense_jobs.append(condense_job)

    # one final job - read the final vcf in, index it, move the index, and write it out
    input_vcf = get_batch().read_input(fragment_files[0])
    final_job = get_batch().new_bash_job(name='index_final_vcf')

    # if there were no jobs, there's no composing...
    if condense_jobs:
        final_job.depends_on(condense_jobs[-1])
    condense_jobs.append(final_job)

    final_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            'vcf.bgz.csi': '{root}.vcf.bgz.csi',
        },
    )
    final_job.image(image_path('bcftools'))
    final_job.storage(config_retrieve(['gcloud_condense', 'storage'], '10Gi'))

    # hypothesis here is that as each separate VCF is block-gzipped, concatenating the results
    # will still be a block-gzipped result. We generate both index types as futureproofing
    final_job.command(f'mv {input_vcf} {final_job.output["vcf.bgz"]}')
    final_job.command(f'tabix {final_job.output["vcf.bgz"]}')
    final_job.command(f'tabix -C {final_job.output["vcf.bgz"]}')

    get_batch().write_output(final_job.output, output_path.removesuffix('.vcf.bgz'))

    # don't start the batch straight away
    return condense_jobs
