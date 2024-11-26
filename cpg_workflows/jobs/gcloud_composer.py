"""
All jobs required to take a shard manifest (containing paths to VCF fragments) and produce a single VCF file
"""

from os.path import join

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, image_path
from cpg_workflows.utils import chunks


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
    chunk_size: int = 30,
) -> tuple[list[str], hb.batch.job.Job]:
    """
    takes a list of things to condense, creates jobs to condense them. Returns a list of the result paths

    Args:
        path_list (list[str]): all the paths to condense
        temp_prefix ():
        chunk_size (): number of objects to condense at once
    Returns:
        list of paths to the condensed objects
    """

    chunk_job = get_batch().new_job(name=f'chunkbuster_{len(path_list)}')
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
    manifest_path: str,
    intermediates_path: str,
    output_path: str,
) -> list[hb.batch.job.Job]:
    """
    compose a series gcloud commands to condense a list of VCF fragments into a single VCF
    gcloud compose has a limit on 32 separate inputs, so at higher input counts we need to do a rolling merge

    compose can only run within a bucket, so we can't do all the intermediate generation in one bucket, and
    then move to a permanent location

    Args:
        manifest_path (str): path to a file in GCP, one GCP file path per line
        intermediates_path (str): path to a bucket, where intermediate files will be written
        output_path (str): path to write the final file to. Must be in the same bucket as the manifest

    Returns:
        a list of the jobs which will generate a single output from the manifest of fragment paths
    """
    manifest_directory = to_path(manifest_path).parent

    # prefix these shard names to get full GCP path for each
    fragment_files = [
        join(manifest_directory, str(fragment).strip()) for fragment in to_path(manifest_path).open().readlines()
    ]

    # rolling squash of the chunks, should enable infinite-ish scaling
    temp_chunk_prefix_num = 1
    condense_jobs: list[hb.batch.job.Job] = []
    while len(fragment_files) > 1:
        condense_temp = join(intermediates_path, f'temp_chunk_{temp_chunk_prefix_num}')
        temp_chunk_prefix_num += 1
        fragment_files, condense_job = compose_condense_fragments(
            fragment_files,
            condense_temp,
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

    final_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
    final_job.image(image_path('bcftools'))
    final_job.storage(config_retrieve(['gcloud_condense', 'storage'], '10Gi'))

    # we specifically need block-gzipped output for some downstream tools
    # otherwise a mv or a bcftools view with --write-index=tbi would be neater
    final_job.command(f'bcftools view {input_vcf} | bgzip -c > {final_job.output["vcf.bgz"]}')
    final_job.command(f'tabix {final_job.output["vcf.bgz"]}')

    get_batch().write_output(final_job.output, output_path.removesuffix('.vcf.bgz'))

    # don't start the batch straight away
    return condense_jobs
