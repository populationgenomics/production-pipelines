#! /usr/bin/env python3


"""
This is a VV simple VEP annotation process, annotating a VCF with a minimal number of annotations
This is to facilitate the migration of the clinvar re-summary process out of AIP/ClinvArbitration
Minimal annotation set to improve runtimes

There are potentially other uses for a VCF -> VCF annotation, but I'm not sure what they'll be
"""
import logging
import os.path
from argparse import ArgumentParser

from hailtop.batch import ResourceFile

from cpg_utils import to_path
from cpg_utils.config import image_path, output_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.bcftools import naive_merge_vcfs

CHROM_LIST: list[str] = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']


def split_vcf_by_chromosome(localised_vcf: ResourceFile, output_dir: str | None = None) -> list[ResourceFile]:
    """
    take a pre-localised file and split it into separate chromosomes
    Args:
        localised_vcf ():
        output_dir (str): if we want these persisted, where to put them

    Returns:

    """
    # split the whole vcf into chromosomes, but keep accurate ordering for all chunks
    ordered_output_vcfs: list[ResourceFile] = []

    # one BCFtools job to rule them all - if no bcftools tasks are required this will still run with no command
    # we shove all subsetting through this one container to reduce VCF copying
    bcftools_job = get_batch().new_job('Subset VCF with bcftools')

    # set some resources
    bcftools_job.image(image_path('bcftools')).cpu(1).memory('8G')

    # todo storage? VCF * 2.5

    # existing fragments, instead of checking for each separately. This is an empty list if the 'dir' doesn't exist yet
    existing_fragments = [str(each_path) for each_path in to_path(output_dir).glob('*')] if output_dir else []
    logging.info(f'Existing VCFs: {existing_fragments}')

    for chrom in CHROM_LIST:

        # if a location to write to was used, check if we already made this (workflow resumption)
        if isinstance(output_dir, str):
            # the name for this chunk of annotation
            result_path = os.path.join(output_dir, f'{chrom}.vcf.bgz')

            # check if it exists already, if so read it in
            if result_path in existing_fragments:
                logging.info(f'{result_path} already exists')
                vcf_fragment = get_batch().read_input_group(
                    **{
                        'vcf.bgz': result_path,
                        'vcf.bgz.tbi': f'{result_path}.tbi',
                    },
                )['vcf.bgz']
                ordered_output_vcfs.append(vcf_fragment)
                continue

        # otherwise declare a new resource group, with a name exclusive to this chromosome
        bcftools_job.declare_resource_group(
            **{chrom: {'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}},
        )

        # create a VCF fragment for this chromosome
        bcftools_job.command(
            f'bcftools view -Oz -o {bcftools_job[chrom]["vcf.bgz"]} -r {chrom} {localised_vcf}',  # type: ignore
        )

        # index it
        bcftools_job.command(f'tabix {bcftools_job[chrom]["vcf.bgz"]} ')  # type: ignore

        if isinstance(output_dir, str):
            # write the fragment & index to GCP
            get_batch().write_output(bcftools_job[chrom], os.path.join(output_dir, chrom))

        # and add to the sorted list for this batch
        ordered_output_vcfs.append(bcftools_job[chrom]['vcf.bgz'])  # type: ignore

    return ordered_output_vcfs


def annotate_localised_vcfs(vcf_list: list[ResourceFile], output_dir: str | None = None) -> list[ResourceFile]:
    """
    annotate each VCF fragment using VEP. Optionally attempt to resume from a folder, and write results to same

    Args:
        vcf_list ():
        output_dir ():

    Returns:
        a list of localised-to-batch annotated VCFs
    """

    # existing outputs
    existing_outputs = [str(each_path) for each_path in to_path(output_dir).glob('*')] if output_dir else []

    logging.info(f'Existing annotated VCFs: {existing_outputs}')

    ordered_annotated: list = []

    # next, annotate!
    for chrom, vcf in zip(CHROM_LIST, vcf_list, strict=True):

        if isinstance(output_dir, str):
            # the name for this chunk of annotation
            result_path = os.path.join(output_dir, f'{chrom}.vcf.bgz')
            # do we already have it generated?
            if result_path in existing_outputs:
                logging.info(f'{result_path} already exists')
                vcf_fragment = get_batch().read_input_group(
                    **{
                        'vcf.bgz': result_path,
                        'vcf.bgz.tbi': f'{result_path}.tbi',
                    },
                )['vcf.bgz']
                ordered_annotated.append(vcf_fragment)
                continue

        # annotate that fragment, making a VCF output
        vep_job = get_batch().new_job(f'annotate {chrom} with VEP')

        # declare a resource group for this annotated VCF output
        vep_job.declare_resource_group(vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # configure the required resources
        vep_job.image(image_path('vep_110')).cpu(4).memory('highmem')

        # gcsfuse works only with the root bucket, without prefix:
        vep_mount_path = to_path(reference_path('vep_110_mount'))
        data_mount = to_path(f'/{vep_mount_path.drive}')
        vep_job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
        vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

        vep_job.command(f'FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz && echo $FASTA')
        vep_job.command(
            f'vep --format vcf --vcf --compress_output bgzip --no_stats --fork 4 --dir_cache {vep_dir}/vep/ '
            f'-o {vep_job.vcf["vcf.bgz"]} '  # type: ignore
            f'-i {vcf} --protein --species homo_sapiens --cache --offline --assembly GRCh38 --fa ${{FASTA}}',
        )
        vep_job.command(f'tabix -p vcf {vep_job.vcf["vcf.bgz"]}')  # type: ignore

        if isinstance(output_dir, str):
            get_batch().write_output(vep_job.vcf, os.path.join(output_dir, chrom))

        # keep a list of the in-batch VCF paths
        ordered_annotated.append(vep_job.vcf["vcf.bgz"])  # type: ignore

    return ordered_annotated


def split_and_annotate_vcf(vcf_in: str | ResourceFile, out_vcf: str):
    """
    take a path to a VCF (in GCP), or a ResourceFile (with implicit accompanying index)
    split that VCF into separate chromosomes (dumb split, but fine for small inputs)
    annotate each chromosomal region separately, then merge the result into a single output

    Args:
        vcf_in (str): path to an input VCF file, or a pre-localised ResourceFile
        out_vcf (str): path to final output
    """

    # allow for situations where this VCF was already a localised Resource
    if isinstance(vcf_in, str):
        vcf_in_batch = get_batch().read_input_group(**{'vcf.bgz': vcf_in, 'vcf.bgz.tbi': f'{vcf_in}.tbi'})['vcf.bgz']
    else:
        vcf_in_batch = vcf_in

    # I want to store the fragmented VCFs in tmp
    vcf_fragments_dir = output_path('vcf_fragments', category='tmp')

    # split the whole vcf into chromosomes, but keep accurate ordering for all chunks
    ordered_output_vcfs: list[ResourceFile] = split_vcf_by_chromosome(
        localised_vcf=vcf_in_batch,
        output_dir=vcf_fragments_dir,
    )

    # new path, also in tmp
    annotated_tmpdir = output_path('annotated_vcf_fragments', category='tmp')

    ordered_annotated: list[ResourceFile] = annotate_localised_vcfs(ordered_output_vcfs, output_dir=annotated_tmpdir)

    # now merge them all
    merged_vcf_in_batch = naive_merge_vcfs(ordered_annotated, output_file=out_vcf, vcfs_localised=True)  # type: ignore
    return merged_vcf_in_batch


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser(description='Run a VCF-in, VCF-out annotation, fragmented by chromosome')
    parser.add_argument('-i', help='VCF in', required=True)
    parser.add_argument('-o', help='VCF out', required=True)
    args = parser.parse_args()

    split_and_annotate_vcf(args.i, out_vcf=args.o)
    get_batch().run(wait=False)
