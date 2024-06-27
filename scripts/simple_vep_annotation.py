#! /usr/bin/env python3


"""
This is a VV simple VEP annotation process, annotating a VCF with a minimal number of annotations
This is to facilitate the migration of the clinvar re-summary process out of AIP/ClinvArbitration
Minimal annotation set to improve runtimes

There are potentially other uses for a VCF -> VCF annotation, but I'm not sure what they'll be
"""
import os.path
from argparse import ArgumentParser

from cpg_utils import to_path
from cpg_utils.config import image_path, output_path, reference_path
from cpg_utils.hail_batch import get_batch


def generate_annotated_data(vcf_in: str):
    """
    if the annotated data Table doesn't exist, generate it

    Args:
        vcf_in (str): path to an input VCF file
    """

    vcf_in_batch = get_batch().read_input_group(**{'vcf.bgz': vcf_in, 'vcf.bgz.tbi': f'{vcf_in}.tbi'})['vcf.bgz']

    # split the whole vcf into chromosomes
    ordered_output_vcfs: list = []

    # One BCFtools job to rule them all
    bcftools_job = get_batch().new_job('Subset VCF with bcftools')

    # set some resources
    bcftools_job.image(image_path('bcftools')).cpu(1).memory('8G')

    # todo storage? VCF * 2.5

    vcf_fragments = output_path('vcf_fragments', category='tmp')

    # existing fragments
    existing_fragments = set(to_path(vcf_fragments).glob('*'))

    for chromosome in [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']:

        # the name for this chunk of annotation
        result_path = os.path.join(vcf_fragments, f'{chromosome}.vcf.bgz')

        # check if it exists
        if result_path in existing_fragments:
            vcf_fragment = get_batch().read_input_group(**{
                'vcf.bgz': result_path,
                'vcf.bgz.tbi': f'{result_path}.tbi'
            })['vcf.bgz']
            ordered_output_vcfs.append(vcf_fragment)
            continue

        # otherwise declare a new resource group, with a name exclusive to this chromosome
        bcftools_job.declare_resource_group(
            **{chromosome: {'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}},
        )

        # create a VCF fragment for this chromosome
        bcftools_job.command(
            f'bcftools view -Oz -o {bcftools_job[chromosome]["vcf.bgz"]} -r {chromosome} {vcf_in_batch}'  # type: ignore
        )
        # index it
        bcftools_job.command(f'tabix {bcftools_job[chromosome]["vcf.bgz"]} ')  # type: ignore

        # write the fragment & index to GCP
        get_batch().write_output(bcftools_job[chromosome], result_path.removesuffix('.vcf.bgz'))

        # and add to the sorted list for this batch
        ordered_output_vcfs.append(bcftools_job[chromosome])

    # new path, not in tmp
    vcf_outputs = output_path('annotated_vcf_fragments')

    # existing outputs
    existing_outputs = set(to_path(vcf_outputs).glob('*'))

    ordered_annotated: list = []

    # next, annotate!
    for chromosome, vcf in zip([f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM'], vcf_fragments):

        # the name for this chunk of annotation
        result_path = os.path.join(vcf_outputs, f'{chromosome}.vcf.bgz')

        if result_path in existing_outputs:
            vcf_fragment = get_batch().read_input_group(**{
                'vcf.bgz': result_path,
                'vcf.bgz.tbi': f'{result_path}.tbi'
            })['vcf.bgz']
            ordered_annotated.append(vcf_fragment)
            continue

        # annotate that fragment, making a VCF output
        vep_job = get_batch().new_job(f'annotate {chromosome} with VEP')

        # declare a resource group for this annotated VCF output
        vep_job.declare_resource_group(vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # configure the required resources
        vep_job.image(image_path('vep_110')).cpu(1).memory('highmem')

        # gcsfuse works only with the root bucket, without prefix:
        vep_mount_path = to_path(reference_path('vep_110_mount'))
        data_mount = to_path(f'/{vep_mount_path.drive}')
        vep_job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
        vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

        vep_job.command(
            f"""\
            FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz
            vep \\
            --format vcf --compress_output bgzip \\
            --vcf \\
            -o {vep_job.vcf} \\ 
            -i {vcf['vcf.bgz']} \\
            --protein \\
            --species homo_sapiens \\ 
            --cache \\
            --offline \\
            --assembly GRCh38 \\
            --dir_cache {vep_dir}/vep/ 
            --fasta $FASTA
            """,
        )

        # this will be an intermediate
        get_batch().write_output(vep_job.vcf, result_path.removesuffix('.vcf.bgz'))

        # # for reasons we need to send the remote paths
        ordered_annotated.append(vep_job.vcf)

    # now merge them all
    # todo

    get_batch().run(wait=False)


if __name__ == '__main__':
    parser = ArgumentParser(description='Run a VCF-in, VCF-out annotation, fragmented by chromosome')
    parser.add_argument('-i', help='VCF in', required=True)
    args = parser.parse_args()

    generate_annotated_data(args.i)
