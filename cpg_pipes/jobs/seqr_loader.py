"""
Jobs specific for seqr-loader.
"""
import logging

from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import Path, images, to_path
from cpg_pipes.hb.batch import hail_query_env
from cpg_pipes.hb.command import wrap_command, python_command
from cpg_pipes.query import seqr_loader
from cpg_pipes.refdata import RefData
from cpg_pipes.types import SequencingType

logger = logging.getLogger(__file__)


def annotate_cohort_jobs(
    b: Batch,
    vcf_path: Path,
    siteonly_vqsr_vcf_path: Path,
    vep_ht_path: Path,
    output_mt_path: Path,
    checkpoints_bucket: Path,
    sequencing_type: SequencingType,
    hail_billing_project: str,
    hail_bucket: Path | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Annotate cohort for seqr loader.
    """
    j = b.new_job(f'annotate cohort', job_attrs)
    j.image(images.DRIVER_IMAGE)
    j.command(
        python_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(vcf_path),
            str(siteonly_vqsr_vcf_path),
            str(vep_ht_path),
            str(output_mt_path),
            overwrite,
            RefData.genome_build,
            sequencing_type.value.upper(),
            str(checkpoints_bucket),
            setup_gcp=True,
            hail_billing_project=hail_billing_project,
            hail_bucket=str(hail_bucket),
            default_reference=RefData.genome_build,
            packages=['cpg_gnomad', 'seqr_loader'],
        )
    )
    return [j]


def annotate_dataset_jobs(
    b: Batch,
    mt_path: Path,
    sample_ids: list[str],
    output_mt_path: Path,
    tmp_bucket: Path,
    hail_billing_project: str,
    hail_bucket: Path | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    subset_mt_path = tmp_bucket / 'cohort-subset.mt'
    subset_j = b.new_job(f'subset cohort to dataset', job_attrs)
    subset_j.image(images.DRIVER_IMAGE)
    subset_j.command(
        python_command(
            seqr_loader,
            seqr_loader.subset_mt_to_samples.__name__,
            str(mt_path),
            sample_ids,
            str(subset_mt_path),
            setup_gcp=True,
            hail_billing_project=hail_billing_project,
            hail_bucket=str(hail_bucket),
            default_reference=RefData.genome_build,
        )
    )

    annotate_j = b.new_job(f'annotate dataset', job_attrs)
    annotate_j.image(images.DRIVER_IMAGE)
    annotate_j.command(
        python_command(
            seqr_loader,
            seqr_loader.annotate_dataset_mt.__name__,
            str(subset_mt_path),
            str(output_mt_path),
            str(tmp_bucket),
            overwrite,
            setup_gcp=True,
            hail_billing_project=hail_billing_project,
            hail_bucket=str(hail_bucket),
            default_reference=RefData.genome_build,
        )
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]


def load_to_es(
    b: Batch,
    mt_path: Path,
    es_host: str,
    es_port: int,
    es_username: str,
    es_password: str,
    es_index: str,
    hail_billing_project: str,
    hail_bucket: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Create a Seqr index for an annotated matrix table.
    """
    # Make a list of dataset samples to subset from the entire matrix table
    j = b.new_job(f'create ES index', job_attrs)
    j.image(images.DRIVER_IMAGE)
    hail_query_env(j, hail_billing_project, hail_bucket)
    cmd = f"""\
    pip3 install click cpg_utils hail seqr_loader elasticsearch
    python3 mt_to_es.py \\
    --mt-path {mt_path} \\
    --es-host {es_host} \\
    --es-port {es_port} \\
    --es-username {es_username} \\
    --es-password {es_password} \\
    --es-index {es_index}
    """
    j.command(
        wrap_command(
            cmd,
            python_script_path=to_path(__file__).parent.parent
            / 'dataproc_scripts'
            / 'seqr'
            / 'mt_to_es.py',
            setup_gcp=True,
        )
    )
    return j
