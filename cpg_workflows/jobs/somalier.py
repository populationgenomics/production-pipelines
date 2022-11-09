"""
Adding jobs for fingerprinting and pedigree checks. Mostly using Somalier.
"""
from typing import cast

import pandas as pd
from hailtop.batch import Batch, Resource
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    image_path,
    reference_path,
    fasta_res_group,
    copy_common_env,
    command,
)
from cpg_workflows.filetypes import CramPath
from cpg_workflows.python_scripts import check_pedigree
from cpg_workflows.resources import STANDARD
from cpg_workflows.targets import Dataset
from cpg_workflows.utils import can_reuse, rich_sample_id_seds

# We want to exclude contaminated samples from relatedness checks. Somalier is not
# designed to work with contaminated samples, and in a presence of contamination it
# can generate a lot of false positive families.
MAX_FREEMIX = 0.04


def pedigree(
    b,
    dataset: Dataset,
    expected_ped_path: Path,
    somalier_path_by_sid: dict[str, Path],
    out_samples_path: Path,
    out_pairs_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    out_checks_path: Path | None = None,
    verifybamid_by_sid: dict[str, Path | str] | None = None,
    label: str | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
) -> list[Job]:
    """
    Add somalier and peddy jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, a path to a fixed PED file if able to recover, and a path to a file
    with relatedness information for each sample pair

    somalier_path_by_sid should be a dict of paths to .somalier fingerprints.
    """
    relate_j = _relate(
        b=b,
        somalier_path_by_sid=somalier_path_by_sid,
        verifybamid_by_sid=verifybamid_by_sid,
        sample_ids=dataset.get_sample_ids(),
        rich_id_map=dataset.rich_id_map(),
        expected_ped_path=expected_ped_path,
        label=label,
        out_samples_path=out_samples_path,
        out_pairs_path=out_pairs_path,
        out_html_path=out_html_path,
        out_html_url=out_html_url,
        job_attrs=job_attrs,
    )

    check_j = _check_pedigree(
        b=b,
        samples_file=relate_j.output_samples,
        pairs_file=relate_j.output_pairs,
        expected_ped=b.read_input(str(expected_ped_path)),
        somalier_html_url=out_html_url,
        rich_id_map=dataset.rich_id_map(),
        dataset_name=dataset.name,
        label=label,
        out_checks_path=out_checks_path,
        job_attrs=job_attrs,
        send_to_slack=send_to_slack,
    )
    check_j.depends_on(relate_j)

    return [relate_j, check_j]


def _make_sample_map(dataset: Dataset):
    """
    Creating sample map to remap internal IDs to participant IDs
    """
    sample_map_fpath = dataset.tmp_prefix() / 'pedigree' / 'sample_map.tsv'
    df = pd.DataFrame(
        [{'id': s.id, 'pid': s.participant_id} for s in dataset.get_samples()]
    )
    if not get_config()['workflow'].get('dry_run', False):
        with sample_map_fpath.open('w') as fp:
            df.to_csv(fp, sep='\t', index=False, header=False)
    return sample_map_fpath


def _check_pedigree(
    b: Batch,
    samples_file: Resource,
    pairs_file: Resource,
    expected_ped: Resource,
    dataset_name: str,
    somalier_html_url: str | None = None,
    rich_id_map: dict[str, str] | None = None,
    label: str | None = None,
    out_checks_path: Path | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
) -> Job:
    """
    Run job that checks pedigree and batch correctness. The job will send a Slack
    message about any mismatches.
    """
    title = 'Pedigree check'
    if label:
        title += f' [{label}]'

    check_j = b.new_job(title, (job_attrs or {}) | dict(tool='python'))
    STANDARD.set_resources(check_j, ncpu=2)
    check_j.image(image_path('cpg_workflows'))

    script_path = to_path(check_pedigree.__file__)
    script_name = script_path.name
    cmd = f"""\
    {rich_sample_id_seds(rich_id_map, [str(samples_file), str(pairs_file), str(expected_ped)])
    if rich_id_map else ''}
    python3 {script_name} \\
    --somalier-samples {samples_file} \\
    --somalier-pairs {pairs_file} \\
    --expected-ped {expected_ped} \\
    {"--html-url {somalier_html_url}" if somalier_html_url else ""} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --{"no-" if not send_to_slack else ""}send-to-slack

    touch {check_j.output}
    """
    if somalier_html_url:
        cmd += f'echo "HTML URL: {somalier_html_url}"'

    copy_common_env(cast(Job, check_j))
    check_j.command(
        command(
            cmd,
            python_script_path=script_path,
            setup_gcp=True,
        )
    )
    if out_checks_path:
        b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _relate(
    b: Batch,
    somalier_path_by_sid: dict[str, Path],
    sample_ids: list[str],
    rich_id_map: dict[str, str],
    expected_ped_path: Path,
    label: str | None,
    out_samples_path: Path,
    out_pairs_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    verifybamid_by_sid: dict[str, Path] | None = None,
    job_attrs: dict | None = None,
) -> Job:
    title = 'Somalier relate'
    if label:
        title += f' [{label}]'

    j = b.new_job(title, (job_attrs or {}) | dict(tool='somalier'))
    j.image(image_path('somalier'))
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    STANDARD.set_resources(j, storage_gb=1 + len(sample_ids) // 4000 * 1)

    cmd = ''
    input_files_file = '$BATCH_TMPDIR/input_files.list'
    samples_ids_file = '$BATCH_TMPDIR/sample_ids.list'
    cmd += f'touch {input_files_file}'
    cmd += f'touch {samples_ids_file}'
    for sample_id in sample_ids:
        if verifybamid_by_sid:
            if sample_id not in verifybamid_by_sid:
                continue
            somalier_file = b.read_input(str(somalier_path_by_sid[sample_id]))
            cmd += f"""
            FREEMIX=$(cat {b.read_input(str(verifybamid_by_sid[sample_id]))} | tail -n1 | cut -f7)
            if [[ $(echo "$FREEMIX > {MAX_FREEMIX}" | bc) -eq 0 ]]; then \
            echo "{somalier_file}" >> {input_files_file}; \
            echo "{sample_id}" >> {samples_ids_file}; \
            fi
            """
        else:
            somalier_file = b.read_input(str(somalier_path_by_sid[sample_id]))
            cmd += f"""
            echo "{somalier_file}" >> {input_files_file}
            echo "{sample_id}" >> {samples_ids_file}
            """

    cmd += f"""
    cat {b.read_input(str(expected_ped_path))} | \
    grep -v Family.ID | \
    grep -f {samples_ids_file} > expected.ped 
    """

    cmd += f"""
    somalier relate \\
    $(cat {input_files_file}) \\
    --ped expected.ped \\
    -o related \\
    --infer
    ls
    mv related.pairs.tsv {j.output_pairs}
    mv related.samples.tsv {j.output_samples}
    {rich_sample_id_seds(rich_id_map, ['related.html'])}
    mv related.html {j.output_html}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'

    j.command(command(cmd))
    # Copy somalier outputs to final destination.
    b.write_output(j.output_samples, str(out_samples_path))
    b.write_output(j.output_pairs, str(out_pairs_path))
    b.write_output(j.output_html, str(out_html_path))
    return j


def extract(
    b,
    cram_path: CramPath,
    out_somalier_path: Path | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = True,
    label: str | None = None,
) -> Job | None:
    """
    Run `somalier extract` to generate a fingerprint (i.e. a `*.somalier` file)
    """
    if can_reuse(out_somalier_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'somalier'}
    j = b.new_job('Somalier extract' + (f' {label}' if label else ''), job_attrs)

    j.image(image_path('somalier'))
    if not cram_path.index_path:
        raise ValueError(
            f'CRAM for somalier is required to have CRAI index ({cram_path})'
        )
    storage_gb = None  # avoid extra disk by default
    if get_config()['workflow']['sequencing_type'] == 'genome':
        storage_gb = 100
    STANDARD.set_resources(j, ncpu=4, storage_gb=storage_gb)

    ref = fasta_res_group(b)
    sites = b.read_input(str(reference_path('somalier_sites')))

    cmd = f"""\
    SITES=$BATCH_TMPDIR/sites/{reference_path('somalier_sites').name}
    retry gsutil cp {reference_path('somalier_sites')} $SITES

    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}
    
    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    somalier extract -d extracted/ --sites {sites} -f {ref.base} $CRAM

    mv extracted/*.somalier {j.output_file}
    """
    j.command(command(cmd, setup_gcp=True, define_retry_function=True))
    b.write_output(j.output_file, str(out_somalier_path))
    return j
