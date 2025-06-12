#! /usr/bin/env python3


"""
This is a simplified VEP annotation process, annotating a VCF with a minimal number of annotations
This is to facilitate the migration of the clinvar re-summary process out of AIP/ClinvArbitration
Minimal annotation set to improve runtimes

Long-running/whole-genome VEP tasks could be improved by localising and unpacking the VEP cache inside the VM
"""
import logging
from argparse import ArgumentParser
from os.path import join

from hailtop.batch import ResourceFile
from hailtop.batch.job import BashJob

from cpg_utils import Path, to_path
from cpg_utils.config import image_path, output_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.bcftools import naive_concat_vcfs

CHROM_LIST: list[str] = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']


VCF_BGZ_SUFFIX: str = 'vcf.bgz'
VCF_BGZ_TBI_SUFFIX: str = 'vcf.bgz.tbi'


def create_vep_job_and_mount_data(job_number: int = 1) -> tuple[BashJob, Path]:
    """

    Args:
        job_number ():

    Returns:

    """
    # annotate that fragment, making a VCF output
    vep_job = get_batch().new_job(f'Annotate part {job_number} with VEP')

    # declare a resource group for this annotated VCF output
    vep_job.declare_resource_group(vcf={VCF_BGZ_SUFFIX: '{root}.vcf.bgz', VCF_BGZ_TBI_SUFFIX: '{root}.vcf.bgz.tbi'})

    # configure the required resources
    vep_job.image(image_path('vep_110')).cpu(4).memory('highmem')

    # gcsfuse works only with the root bucket, without prefix:
    vep_mount_path = to_path(reference_path('vep_110_mount'))
    data_mount = to_path(f'/{vep_mount_path.drive}')
    vep_job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])
    return vep_job, vep_dir


def split_vcf_by_chromosome(
    localised_vcf: ResourceFile,
    output_dir: str | None = None,
    storage: str | None = None,
) -> tuple[list[ResourceFile], BashJob]:
    """
    take a pre-localised file and split it into separate chromosomes
    Args:
        localised_vcf ():
        output_dir (str): if we want these persisted, where to put them
        storage (str): can be set if required, default OS partition can hold ~5GB without additional attached storage

    Returns:
        the list of resource files, and the job that created them
    """
    # split the whole vcf into chromosomes, but keep accurate ordering for all chunks
    ordered_output_vcfs: list[ResourceFile] = []

    # one BCFtools job to rule them all - if no bcftools tasks are required this will still run with no command
    # we shove all subsetting through this one container to reduce VCF copying
    bcftools_job = get_batch().new_bash_job('Subset VCF with bcftools')

    # set some resources
    bcftools_job.image(image_path('bcftools')).cpu(1).memory('8G')

    # only attach further storage if requested
    if storage:
        bcftools_job.storage(storage)

    # existing fragments, instead of checking for each separately. This is an empty list if the 'dir' doesn't exist yet
    existing_fragments = [str(each_path) for each_path in to_path(output_dir).glob('*')] if output_dir else []
    logging.info(f'Existing VCFs: {existing_fragments}')

    for chrom in CHROM_LIST:

        # if a location to write to was used, check if we already made this (workflow resumption)
        if isinstance(output_dir, str):
            # the name for this chunk of annotation
            result_path = join(output_dir, f'{chrom}.vcf.bgz')

            # check if it exists already, if so read it in
            if result_path in existing_fragments:
                logging.info(f'{result_path} already exists')
                vcf_fragment = get_batch().read_input_group(
                    **{
                        VCF_BGZ_SUFFIX: result_path,
                        VCF_BGZ_TBI_SUFFIX: f'{result_path}.tbi',
                    },
                )[VCF_BGZ_SUFFIX]
                ordered_output_vcfs.append(vcf_fragment)
                continue

        # otherwise declare a new resource group, with a name exclusive to this chromosome
        bcftools_job.declare_resource_group(
            **{chrom: {VCF_BGZ_SUFFIX: '{root}.vcf.bgz', VCF_BGZ_TBI_SUFFIX: '{root}.vcf.bgz.tbi'}},
        )

        # create a VCF fragment for this chromosome, and index the result
        bcftools_job.command(
            f'bcftools view -Oz --write-index=tbi -o {bcftools_job[chrom][VCF_BGZ_SUFFIX]} -r {chrom} {localised_vcf}',
        )

        if isinstance(output_dir, str):
            # write the fragment & index to GCP
            get_batch().write_output(bcftools_job[chrom], join(output_dir, chrom))

        # and add to the sorted list for this batch
        ordered_output_vcfs.append(bcftools_job[chrom][VCF_BGZ_SUFFIX])

    return ordered_output_vcfs, bcftools_job


def minimal_annotation(
    vcf_list: list[ResourceFile],
    output_dir: str,
) -> tuple[list[ResourceFile], list[BashJob]]:
    """
    annotate each VCF fragment using VEP. Optionally attempt to resume from a folder, and write results to same

    I've stashed a tarball of the VEP 110 cache here: gs://cpg-common-test/references/vep_110_compressed.tar.gz
    instead of enabling a cloud-fused annotation source, we could copy this into the VM, unpack, and go...

    vep_110_tarball = get_batch().read_input('gs://cpg-common-test/references/vep_110_compressed.tar.gz')
    job.command('tar -xzf compressed_vep.tar.gz ${{BATCH_TMPDIR}}/vep')
    --dir_cache ${{BATCH_TMPDIR}}/vep --fa ${{BATCH_TMPDIR}}/vep/homo_sapiens/110/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    For the clinvar use-case the runtime isn't enough to merit this optimisation, but for a genome-wide annotation in
    a single VM, probably best to localise the files?

    Args:
        vcf_list ():
        output_dir ():

    Returns:
        a list of localised-to-batch annotated VCFs, and the annotation jobs
    """
    existing_outputs = [str(each_path) for each_path in to_path(output_dir).glob('*')] if output_dir else []

    logging.info(f'Existing annotated VCFs: {existing_outputs}')

    ordered_annotated: list = []
    annotation_jobs: list = []

    # next, annotate!
    for job_number, vcf in enumerate(vcf_list, start=1):

        # the name for this chunk of annotation
        result_path = join(output_dir, f'{job_number}.vcf.bgz')
        # do we already have it generated?
        if result_path in existing_outputs:
            logging.info(f'{result_path} already exists')
            vcf_fragment = get_batch().read_input_group(
                **{VCF_BGZ_SUFFIX: result_path, VCF_BGZ_TBI_SUFFIX: f'{result_path}.tbi'},
            )[VCF_BGZ_SUFFIX]
            ordered_annotated.append(vcf_fragment)
            continue

        # annotate that fragment, making a VCF output
        vep_job, vep_dir = create_vep_job_and_mount_data(job_number)

        vep_job.command(f'FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz && echo $FASTA')
        vep_job.command(
            f"""
            vep \
            -i {vcf} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --no_stats \
            --fork 4 \
            --dir_cache {vep_dir}/vep/ \
            -o {vep_job.vcf[VCF_BGZ_SUFFIX]} \
            --protein \
            --species homo_sapiens \
            --cache \
            --offline \
            --assembly GRCh38 \
            --fa ${{FASTA}}
            """,
        )
        vep_job.command(f'tabix -p vcf {vep_job.vcf[VCF_BGZ_SUFFIX]}')
        get_batch().write_output(vep_job.vcf, join(output_dir, str(job_number)))

        annotation_jobs.append(vep_job)
        # keep a list of the in-batch VCF paths
        ordered_annotated.append(vep_job.vcf[VCF_BGZ_SUFFIX])

    return ordered_annotated, annotation_jobs


def annotate_localised_vcfs(
    vcf_list: list[ResourceFile],
    output_dir: str,
) -> tuple[list[ResourceFile], list[BashJob]]:
    """
    a slightly beefier method for annotating VCFs, based on the minimal method above
    this is used to apply the general annotations used in Talos, and can be used as a
    template for VEP installation and annotation

    Args:
        vcf_list ():
        output_dir ():

    Returns:
        a list of localised-to-batch annotated VCFs, and the annotation jobs
    """

    # existing outputs
    existing_outputs = [str(each_path) for each_path in to_path(output_dir).glob('*')] if output_dir else []

    logging.info(f'Existing annotated VCFs: {existing_outputs}')

    ordered_annotated: list = []
    annotation_jobs: list = []

    # next, annotate!
    for job_number, vcf in enumerate(vcf_list, start=1):

        # the name for this chunk of annotation
        result_path = join(output_dir, f'{job_number}.vcf.bgz')
        # do we already have it generated?
        if result_path in existing_outputs:
            logging.info(f'{result_path} already exists')
            vcf_fragment = get_batch().read_input_group(
                **{VCF_BGZ_SUFFIX: result_path, VCF_BGZ_TBI_SUFFIX: f'{result_path}.tbi'},
            )[VCF_BGZ_SUFFIX]
            ordered_annotated.append(vcf_fragment)
            continue

        # annotate that fragment, making a VCF output
        vep_job, vep_dir = create_vep_job_and_mount_data(job_number)

        loftee_conf = {
            'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
            'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
            'conservation_file': f'{vep_dir}/loftee.sql',
            'loftee_path': '$VEP_DIR_PLUGINS',
        }

        vep_job.command(f'FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz && echo $FASTA')
        vep_job.command(
            f"""
            vep \
            --format vcf \
            -i {vcf} \
            --vcf \
            --compress_output bgzip \
            --no_stats \
            --dir_cache {vep_dir}/vep/ \
            --species homo_sapiens \
            --cache \
            --offline \
            --assembly GRCh38 \
            --fa ${{FASTA}} \
            -o {vep_job.vcf[VCF_BGZ_SUFFIX]} \
            --protein \
            --af_gnomadg \
            --af_gnomade \
            --mane_select \
            --plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz \
            --plugin LoF,{",".join(f"{k}:{v}" for k, v in loftee_conf.items())} \
            --plugin UTRAnnotator,file=$UTR38
            """,
        )
        vep_job.command(f'tabix -p vcf {vep_job.vcf[VCF_BGZ_SUFFIX]}')
        get_batch().write_output(vep_job.vcf, join(output_dir, str(job_number)))

        annotation_jobs.append(vep_job)
        # keep a list of the in-batch VCF paths
        ordered_annotated.append(vep_job.vcf[VCF_BGZ_SUFFIX])

    return ordered_annotated, annotation_jobs


def split_and_annotate_vcf(
    vcf_in: str | ResourceFile,
    out_vcf: str,
    minimal: bool = True,
) -> tuple[ResourceFile, list[BashJob]]:
    """
    take a path to a VCF (in GCP), or a ResourceFile (with implicit accompanying index)
    split that VCF into separate chromosomes (dumb split, but fine for small inputs)
    annotate each chromosomal region separately, then merge the result into a single output

    Args:
        vcf_in (str): path to an input VCF file, or a pre-localised ResourceFile
        out_vcf (str): path to final output
        minimal (bool): whether to use the minimal annotation set, default is True (for ClinvArbitration)
    """

    # allow for situations where this VCF was already a localised Resource
    if isinstance(vcf_in, str):
        vcf_in_batch = get_batch().read_input_group(**{VCF_BGZ_SUFFIX: vcf_in, VCF_BGZ_TBI_SUFFIX: f'{vcf_in}.tbi'})[
            VCF_BGZ_SUFFIX
        ]
    else:
        vcf_in_batch = vcf_in

    # I want to store the fragmented VCFs in tmp
    vcf_fragments_dir = output_path('vcf_fragments', category='tmp')

    # split the whole vcf into chromosomes, but keep accurate ordering for all chunks
    ordered_output_vcfs, bcftools_job = split_vcf_by_chromosome(
        localised_vcf=vcf_in_batch,
        output_dir=vcf_fragments_dir,
    )

    # new path, also in tmp
    annotated_tmpdir = output_path('annotated_vcf_fragments', category='tmp')

    # choose the annotation method
    annotation_method = minimal_annotation if minimal else annotate_localised_vcfs

    ordered_annotated, annotation_jobs = annotation_method(ordered_output_vcfs, output_dir=annotated_tmpdir)

    # now merge them all
    merged_vcf_in_batch, merge_job = naive_concat_vcfs(ordered_annotated, output_file=out_vcf, vcfs_localised=True)

    # collect all the jobs
    all_jobs = [bcftools_job, *annotation_jobs, merge_job]
    return merged_vcf_in_batch, all_jobs


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser(description='Run a VCF-in, VCF-out annotation, fragmented by chromosome')
    parser.add_argument('-i', help='VCF in', required=True)
    parser.add_argument('-o', help='VCF out', required=True)
    parser.add_argument('-m', help='Use to request most minimal annotation', action='store_true')
    args = parser.parse_args()

    _resource_file, _all_jobs = split_and_annotate_vcf(args.i, out_vcf=args.o, minimal=args.m)
    get_batch().run(wait=False)
