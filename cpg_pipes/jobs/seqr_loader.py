"""
Jobs specific for seqr-loader.
"""
from hailtop.batch.job import Job

from cpg_pipes import Path, images, utils
from cpg_pipes.hb.batch import Batch, hail_query_env
from cpg_pipes.hb.command import wrap_command


def annotate_dataset(
    b: Batch,
    annotated_mt_path: Path,
    sample_ids: list[str],
    output_mt_path: Path,
    tmp_bucket: Path,
    hail_billing_project: str,
    hail_bucket: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr)
    """
    # Make a list of dataset samples to subset from the entire matrix table
    subset_path = tmp_bucket / 'seqr-samples.txt'
    with subset_path.open('w') as f:
        f.write('\n'.join(sample_ids))

    j = b.new_job(f'split and annotate', job_attrs)
    j.image(images.DRIVER_IMAGE)
    hail_query_env(j, hail_billing_project, hail_bucket)
    cmd = f"""\
    pip3 install click cpg_utils hail seqr_loader
    python3 subset_mt.py \\
    --mt-path {annotated_mt_path} \\
    --out-mt-path {output_mt_path} \\
    --subset-tsv {subset_path}
    """
    j.command(wrap_command(
        cmd,
        python_script=utils.QUERY_SCRIPTS_DIR / 'seqr' / 'subset_mt.py',
        setup_gcp=True,
    ))
    return j


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
    pip3 install click cpg_utils hail seqr_loader elasticsearch==7.9.1
    python3 mt_to_es.py \\
    --mt-path {mt_path} \\
    --es-host {es_host} \\
    --es-port {es_port} \\
    --es-username {es_username} \\
    --es-password {es_password} \\
    --es-index {es_index}
    """
    j.command(wrap_command(
        cmd,
        python_script=utils.QUERY_SCRIPTS_DIR / 'seqr' / 'mt_to_es.py',
        setup_gcp=True,
    ))
    return j
