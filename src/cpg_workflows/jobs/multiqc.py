#!/usr/bin/env python3

"""
Batch jobs to run MultiQC.
"""
from typing import cast

from hailtop.batch import Batch, ResourceFile
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import command, copy_common_env
from cpg_workflows.python_scripts import check_multiqc
from cpg_workflows.resources import STANDARD
from cpg_workflows.targets import Dataset
from cpg_workflows.utils import rich_sequencing_group_id_seds


def multiqc(
    b: Batch,
    dataset: Dataset,
    tmp_prefix: Path,
    paths: list[Path],
    out_json_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    out_checks_path: Path | None = None,
    label: str | None = None,
    ending_to_trim: set[str] | None = None,
    modules_to_trim_endings: set[str] | None = None,
    job_attrs: dict | None = None,
    sequencing_group_id_map: dict[str, str] | None = None,
    extra_config: dict | None = None,
) -> list[Job]:
    """
    Run MultiQC for the files in `qc_paths`
    @param b: batch object
    @param tmp_prefix: bucket for tmp files
    @param paths: file bucket paths to pass into MultiQC
    @param dataset: Dataset object
    @param out_json_path: where to write MultiQC-generated JSON file
    @param out_html_path: where to write the HTML report
    @param out_html_url: URL corresponding to the HTML report
    @param out_checks_path: flag indicating that QC checks were done
    @param label: To add to the report's, Batch job's, and Slack message's titles
    @param ending_to_trim: trim these endings from input files to get sequencing group names
    @param modules_to_trim_endings: list of modules for which trim the endings
    @param job_attrs: attributes to add to Hail Batch job
    @param sequencing_group_id_map: sequencing group ID map for bulk sequencing group renaming:
        (https://multiqc.info/docs/#bulk-sample-renaming-in-reports)
    @param extra_config: extra config to pass to MultiQC
    @return: job objects
    """
    title = 'MultiQC'
    if label:
        title += f' [{label}]'

    mqc_j = b.new_job(title, (job_attrs or {}) | dict(tool='MultiQC'))
    mqc_j.image(image_path('multiqc'))
    STANDARD.set_resources(mqc_j, ncpu=16)

    file_list_path = tmp_prefix / f'{dataset.alignment_inputs_hash()}_multiqc-file-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in paths])
    file_list = b.read_input(str(file_list_path))

    endings_conf = ', '.join(list(ending_to_trim)) if ending_to_trim else ''
    modules_conf = ', '.join(list(modules_to_trim_endings)) if modules_to_trim_endings else ''

    if sequencing_group_id_map:
        sample_map_path = tmp_prefix / f'{dataset.alignment_inputs_hash()}_rename-sample-map.tsv'
        if not get_config()['workflow'].get('dry_run', False):
            _write_sg_id_map(sequencing_group_id_map, sample_map_path)
        sample_map_file = b.read_input(str(sample_map_path))
    else:
        sample_map_file = None

    if extra_config:
        serialised = ', '.join(f'{k}: {v}' for k, v in extra_config.items())
        extra_config_param = f'--cl-config "{serialised}"'
    else:
        extra_config_param = ''

    report_filename = 'report'
    cmd = f"""\
    mkdir inputs
    cat {file_list} | gsutil -m cp -I inputs/

    multiqc -f inputs -o output \\
    {f"--replace-names {sample_map_file} " if sample_map_file else ''} \\
    --title "{title} for dataset <b>{dataset.name}</b>" \\
    --filename {report_filename}.html \\
    --cl-config "extra_fn_clean_exts: [{endings_conf}]" \\
    --cl-config "max_table_rows: 10000" \\
    --cl-config "use_filename_as_sample_name: [{modules_conf}]" \\
    {extra_config_param}

    ls output/{report_filename}_data
    cp output/{report_filename}.html {mqc_j.html}
    cp output/{report_filename}_data/multiqc_data.json {mqc_j.json}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'

    mqc_j.command(command(cmd, setup_gcp=True))
    b.write_output(mqc_j.html, str(out_html_path))
    b.write_output(mqc_j.json, str(out_json_path))

    assert isinstance(mqc_j.json, ResourceFile)
    jobs: list[Job] = [mqc_j]
    if get_config().get('qc_thresholds'):
        check_j = check_report_job(
            b=b,
            multiqc_json_file=mqc_j.json,
            multiqc_html_url=out_html_url,
            rich_id_map=dataset.rich_id_map(),
            dataset_name=dataset.name,
            label=label,
            out_checks_path=out_checks_path,
            job_attrs=job_attrs,
        )
        check_j.depends_on(mqc_j)
        jobs.append(check_j)
    return jobs


def check_report_job(
    b: Batch,
    multiqc_json_file: ResourceFile,
    dataset_name: str,
    multiqc_html_url: str | None = None,
    label: str | None = None,
    rich_id_map: dict[str, str] | None = None,
    out_checks_path: Path | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
) -> Job:
    """
    Run job that checks MultiQC JSON result and sends a Slack notification
    about failed samples.
    """
    title = 'MultiQC'
    if label:
        title += f' [{label}]'
    check_j = b.new_job(f'{title} check', (job_attrs or {}) | dict(tool='python'))
    STANDARD.set_resources(check_j, ncpu=2)
    check_j.image(image_path('cpg_workflows'))

    script_path = to_path(check_multiqc.__file__)
    script_name = script_path.name
    cmd = f"""\
    {rich_sequencing_group_id_seds(rich_id_map, [multiqc_json_file])
    if rich_id_map else ''}

    python3 {script_name} \\
    --multiqc-json {multiqc_json_file} \\
    --html-url {multiqc_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --{"no-" if not send_to_slack else ""}send-to-slack

    touch {check_j.output}
    echo "HTML URL: {multiqc_html_url}"
    """

    copy_common_env(cast(Job, check_j))
    check_j.command(
        command(
            cmd,
            python_script_path=script_path,
            setup_gcp=True,
        ),
    )
    if out_checks_path:
        b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _write_sg_id_map(sequencing_group_map: dict[str, str], out_path: Path):
    """
    Configuring MultiQC to support bulk sequencing group rename. `sequencing_group_map` is a dictionary
    of sequencing group IDs. The map doesn't have to have records for all sequencing groups.
    Example:
    {
        'SID1': 'Patient1',
        'SID2': 'Patient2'
    }
    https://multiqc.info/docs/#bulk-sample-renaming-in-reports
    """
    with out_path.open('w') as fh:
        for sgid, new_sgid in sequencing_group_map.items():
            fh.write('\t'.join([sgid, new_sgid]) + '\n')
