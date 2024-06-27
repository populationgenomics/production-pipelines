"""
general tasks using bcftools in a pipeline context
"""

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch


def naive_merge_vcfs(
    input_list: list[str],
    output_file: str,
    cpu: int = 4,
    memory: str = '16Gi',
    storage: str = '50Gi',
    missing_to_ref: bool = False,
    vcfs_localised: bool = False
):
    """
    a generic method for taking multiple vcfs, merging into a single vcf

    Args:
        input_list (list[str]): all VCFs to merge, these will be read into the batch and merged
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
        f'{merge_job.output["vcf.bgz"]} --threads {cpu} -m all {" -0" if missing_to_ref else ""}',  # type: ignore
    )
    merge_job.command(f'tabix {merge_job.output["vcf.bgz"]}')  # type: ignore

    # write the result out
    get_batch().write_output(merge_job.output, output_file.removesuffix('vcf.bgz'))

    return merge_job.output["vcf.bgz"]  # type: ignore
