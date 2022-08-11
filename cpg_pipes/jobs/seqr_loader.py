"""
Jobs specific for seqr-loader.
"""
import logging

from cpg_utils.hail_batch import image_path, genome_build, reference_path
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import Path, to_path
from cpg_pipes.hb.command import wrap_command, python_command
from cpg_pipes.query import seqr_loader

logger = logging.getLogger(__file__)


def annotate_cohort_jobs(
    b: Batch,
    vcf_path: Path,
    siteonly_vqsr_vcf_path: Path,
    vep_ht_path: Path,
    output_mt_path: Path,
    checkpoint_prefix: Path,
    sequencing_type: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Annotate cohort for seqr loader.
    """
    j = b.new_job(f'annotate cohort', job_attrs)
    j.image(image_path('hail'))
    j.command(
        python_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(vcf_path),
            str(siteonly_vqsr_vcf_path),
            str(vep_ht_path),
            str(output_mt_path),
            overwrite,
            genome_build(),
            sequencing_type,
            str(checkpoint_prefix),
            setup_gcp=True,
            setup_hail=True,
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
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    subset_mt_path = tmp_bucket / 'cohort-subset.mt'
    subset_j = b.new_job(
        f'subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'}
    )
    subset_j.image(image_path('hail'))
    assert sample_ids
    subset_j.command(
        python_command(
            seqr_loader,
            seqr_loader.subset_mt_to_samples.__name__,
            str(mt_path),
            sample_ids,
            str(subset_mt_path),
            setup_gcp=True,
            setup_hail=True,
        )
    )

    annotate_j = b.new_job(
        f'annotate dataset', (job_attrs or {}) | {'tool': 'hail query'}
    )
    annotate_j.image(image_path('hail'))
    annotate_j.command(
        python_command(
            seqr_loader,
            seqr_loader.annotate_dataset_mt.__name__,
            str(subset_mt_path),
            str(output_mt_path),
            str(tmp_bucket),
            overwrite,
            setup_gcp=True,
            setup_hail=True,
        )
    )
    annotate_j.depends_on(subset_j)
    return [subset_j, annotate_j]


def load_to_es(
    b: Batch,
    mt_path: Path,
    es_index: str,
    job_attrs: dict | None = None,
) -> Job:
    """
    Create a Seqr index for an annotated matrix table.
    """
    # Make a list of dataset samples to subset from the entire matrix table
    j = b.new_job(f'create ES index', (job_attrs or {}) | {'tool': 'hail query'})
    j.image(image_path('hail'))
    cmd = f"""\
    pip3 install click cpg_utils hail seqr_loader elasticsearch
    python3 mt_to_es.py \\
    --mt-path {mt_path} \\
    --es-index {es_index} \\
    --liftover-path {reference_path('liftover_38_to_37')}
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
