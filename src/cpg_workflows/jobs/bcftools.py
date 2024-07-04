"""
general tasks using bcftools in a pipeline context
"""

from hailtop.batch import ResourceFile
from hailtop.batch.job import BashJob

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

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
        f'{merge_job.output["vcf.bgz"]} --threads {cpu} -m all {" -0" if missing_to_ref else ""} --write-index=tbi',
    )

    # write the result out
    get_batch().write_output(merge_job.output, output_file.removesuffix('.vcf.bgz'))

    return merge_job.output["vcf.bgz"], merge_job
