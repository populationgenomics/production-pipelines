"""
jobs relating to the validation steps of the pipeline
"""


from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, query_command


def validation_mt_to_vcf_job(
    b: Batch,
    mt_path: Path,
    sample_id: str,
    out_vcf_path: Path,
    job_attrs: dict | None = None,
    depends_on: list[Job] | None = None
):
    """
    Take the single-dataset MT, and write to a VCF

    Args:
        b (hb.Batch): the batch to add jobs into
        mt_path (str): path of the AnnotateDataset MT
        sample_id (str): sample name
        out_vcf_path (str): path to write new VCF to
        job_attrs (dict):
        depends_on (hb.Job|list[hb.Job]): jobs to depend on

    Returns:
        this single Job
    """
    from cpg_workflows.query_modules import validation

    vcf_j = b.new_job(
        f'VCF from dataset MT', (job_attrs or {}) | {'tool': 'hail query'}
    )
    vcf_j.image(image_path('cpg_workflows'))
    vcf_j.command(
        query_command(
            validation,
            validation.single_sample_vcf_from_dataset_vcf.__name__,
            str(mt_path),
            sample_id,
            out_vcf_path,
            setup_gcp=True,
        )
    )
    if depends_on:
        vcf_j.depends_on(*depends_on)

    return vcf_j
