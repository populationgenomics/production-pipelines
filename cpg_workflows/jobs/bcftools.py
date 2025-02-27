"""
general tasks using bcftools in a pipeline context
"""

from hailtop.batch import ResourceFile
from hailtop.batch.job import BashJob

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import GvcfPath

# mypy: ignore_errors


def naive_concat_vcfs(
    input_list: list[str | ResourceFile],
    output_file: str,
    cpu: int = 4,
    memory: str = '16Gi',
    storage: str = '50Gi',
    vcfs_localised: bool = False,
) -> tuple[ResourceFile, BashJob]:
    """
    A generic method to concatenate multiple VCFs with identical sample sets to create one VCF.
    Sample IDs must be identical in each component VCF, or all missing (sites-only)

    Args:
        input_list (list[str]): all VCFs to concatenate, can be pre-localised (see vcfs_localised[bool)
        output_file (str): path to a vcf.bgz file to write to
        cpu (int): number of cores (threads when merging)
        memory (str): RAM requirement
        storage (str): storage requirement for the task
        vcfs_localised (bool): if false, read into batch. If true, assume already in the batch

    Returns:
        the ResourceFile in-batch for the merged VCF
    """
    if vcfs_localised:
        batch_vcfs = input_list
    else:
        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in input_list
        ]

    concat = get_batch().new_job('Concat VCFs', attributes={'tool': 'bcftools'})
    concat.image(image=image_path('bcftools'))

    # guessing at resource requirements
    concat.cpu(cpu)
    concat.memory(memory)
    concat.storage(storage)
    concat.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # --write-index: create tabix index
    # -n: concatenate without recompression (not used, cannot be combined with --write-index)
    concat.command(
        f'bcftools concat {" ".join(batch_vcfs)} -Oz -o {concat.output["vcf.bgz"]} --threads {cpu} --write-index=tbi',
    )

    # write the result out
    get_batch().write_output(concat.output, output_file.removesuffix('.vcf.bgz'))

    return concat.output["vcf.bgz"], concat


def naive_merge_vcfs(
    input_list: list[str | ResourceFile],
    output_file: str,
    cpu: int = 4,
    memory: str = '16Gi',
    storage: str = '50Gi',
    missing_to_ref: bool = False,
    vcfs_localised: bool = False,
) -> tuple[ResourceFile, BashJob]:
    """
    A generic method to merge multiple VCFs with non-overlapping sample sets to create one multi-sample VCF.
    Sample names should be unique across all files.

    Args:
        input_list (list[str]): all VCFs to merge, can be pre-localised (see vcfs_localised[bool)
        output_file (str): path to a vcf.bgz file to write to
        cpu (int): number of cores (threads when merging)
        memory (str): RAM requirement
        storage (str): storage requirement for the task
        missing_to_ref (bool): if specified, replace missing calls with unphased WT 0/0
        vcfs_localised (bool): if false, read into batch. If true, assume already in the batch

    Returns:
        the ResourceFile in-batch for the merged VCF
    """

    if vcfs_localised:
        batch_vcfs = input_list
    else:
        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in input_list
        ]

    merge_job = get_batch().new_job('Merge VCFs', attributes={'tool': 'bcftools'})
    merge_job.image(image=image_path('bcftools'))

    # guessing at resource requirements
    merge_job.cpu(cpu)
    merge_job.memory(memory)
    merge_job.storage(storage)
    merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: missing-calls-to-ref (not used by default)
    merge_job.command(
        f'bcftools merge {" ".join(batch_vcfs)} -Oz -o '
        f'{merge_job.output["vcf.bgz"]} --threads {cpu} -m none {" -0" if missing_to_ref else ""} --write-index=tbi',
    )

    # write the result out
    get_batch().write_output(merge_job.output, output_file.removesuffix('.vcf.bgz'))

    return merge_job.output["vcf.bgz"], merge_job


def merge_ss_vcfs(
    input_list: list[str],
    output_file: str,
    cpu: int = 4,
    memory: str = '16Gi',
    storage: str = '50Gi',
    missing_to_ref: bool = False,
    region_file: str | None = None,
) -> BashJob:
    """
    A generic method to merge multiple VCFs with non-overlapping sample sets to create one multi-sample VCF.
    Sample names should be unique across all files.

    Args:
        input_list (list[str]): all VCFs to merge, can be pre-localised (see vcfs_localised[bool)
        output_file (str): path to a vcf.bgz file to write to
        cpu (int): number of cores (threads when merging)
        memory (str): RAM requirement
        storage (str): storage requirement for the task
        missing_to_ref (bool): if specified, replace missing calls with unphased WT 0/0
        region_file (str): path to a BED file with the region to extract

    Returns:
        the ResourceFile in-batch for the merged VCF
    """

    batch_vcfs = [
        get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
        for each_vcf in input_list
    ]

    job = get_batch().new_job('Merge SS VCFs', attributes={'tool': 'bcftools'})
    job.image(image=image_path('bcftools_120'))

    # guessing at resource requirements
    job.cpu(cpu)
    job.memory(memory)
    job.storage(storage)
    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    region_string = ''
    if region_file:
        localised_region_file = get_batch().read_input(region_file)
        region_string = f' -R {localised_region_file}'

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy (don't create multiallelic variants)
    # -0: missing-calls-to-ref (not used by default)
    # -R: region to extract (optional depending on the region file being passed as an argument)
    job.command(
        f'bcftools merge {" ".join(batch_vcfs)} {region_string} -Oz -o temp.vcf.gz '
        f'--threads {cpu} -m none {" -0" if missing_to_ref else ""} --write-index=tbi ',
    )
    # use the fill-tags plugin to correct the combined callset AC/AN
    job.command(f'bcftools +fill-tags temp.vcf.gz -Oz -o {job.output["vcf.bgz"]} --write-index=tbi -- -t AF ')

    # write the result out
    get_batch().write_output(job.output, output_file.removesuffix('.vcf.bgz'))

    return job


def strip_gvcf_to_vcf(
    gvcf: GvcfPath,
    output: str,
    job_attrs: dict | None = None,
) -> BashJob:
    """
    Takes a GVCF and writes out a VCF
    Args:
        gvcf (GvcfPath): the GVCF to strip
        output (str): path to write the VCF
        job_attrs (dict): job attributes
    """

    # read in the gVCF (no need for an index)
    localised_gvcf = get_batch().read_input(str(gvcf.path))

    job = get_batch().new_job('Strip GVCF to VCF', job_attrs)

    # declare a resource group for the output
    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # use a bcftools which can write an inline index
    job.image(image_path('bcftools_120'))

    job.command(
        f'bcftools view -m3 {localised_gvcf} | '
        f'bcftools norm -m -any | '
        f'grep -v NON_REF | '
        f'bcftools view -Oz -o {job.output["vcf.bgz"]} --write-index=tbi',
    )
    get_batch().write_output(job.output, output.removesuffix('.vcf.bgz'))
    return job
