#!/usr/bin/env python3

"""
Batch pipeline to run WGS QC.
"""

import logging

from cpg_utils.hail_batch import image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import Path
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.slack import slack_message_cmd

logger = logging.getLogger(__file__)


def multiqc(
    b: Batch,
    tmp_prefix: Path,
    paths: list[Path],
    dataset_name: str,
    out_html_path: Path,
    out_json_path: Path,
    out_html_url: str | None = None,
    ending_to_trim: set[str] | None = None,
    modules_to_trim_endings: set[str] | None = None,
    job_attrs: dict | None = None,
    sample_id_map: dict[str, str] | None = None,
    extra_config: dict | None = None,
) -> Job:
    """
    Run MultiQC for the files in `qc_paths`
    @param b: batch object
    @param tmp_prefix: bucket for tmp files
    @param paths: file bucket paths to pass into MultiQC
    @param dataset_name: dataset name
    @param out_json_path: where to write MultiQC-generated JSON file
    @param out_html_path: where to write the HTML report
    @param out_html_url: URL corresponding to the HTML report
    @param ending_to_trim: trim these endings from input files to get sample names
    @param modules_to_trim_endings: list of modules for which trim the endings
    @param job_attrs: attributes to add to Hail Batch job
    @param sample_id_map: sample ID map for bulk sample renaming:
        (https://multiqc.info/docs/#bulk-sample-renaming-in-reports)
    @param extra_config: extra config to pass to MultiQC
    @return: job object
    """
    j = b.new_job('Run MultiQC', job_attrs)
    j.image(image_path('multiqc'))
    STANDARD.set_resources(j, ncpu=16)

    file_list_path = tmp_prefix / 'multiqc-file-list.txt'
    with file_list_path.open('w') as f:
        f.writelines([f'{p}\n' for p in paths])
    file_list = b.read_input(str(file_list_path))

    endings_conf = ', '.join(list(ending_to_trim)) if ending_to_trim else ''
    modules_conf = (
        ', '.join(list(modules_to_trim_endings)) if modules_to_trim_endings else ''
    )

    if sample_id_map:
        sample_map_path = tmp_prefix / 'rename-sample-map.tsv'
        _write_sample_id_map(sample_id_map, sample_map_path)
        sample_map_file = b.read_input(str(sample_map_path))
    else:
        sample_map_file = None

    if extra_config:
        serialised = ", ".join(f"{k}: {v}" for k, v in extra_config.items())
        extra_config_param = f'--cl-config "{serialised}"'
    else:
        extra_config_param = ''

    report_filename = 'report'
    cmd = f"""\
    mkdir inputs
    cat {file_list} | gsutil -m cp -I inputs/
    
    # Temporary fix for Somalier module before https://github.com/ewels/MultiQC/pull/1670
    # is merged and released:
    git clone --single-branch --branch fix-somalier https://github.com/vladsaveliev/MultiQC
    pip3 install -e MultiQC

    multiqc -f inputs -o output \\
    {f"--replace-names {sample_map_file} " if sample_map_file else ''} \\
    --title "CRAM QC report for dataset {dataset_name}" \\
    --filename {report_filename}.html \\
    --cl-config "extra_fn_clean_exts: [{endings_conf}]" \\
    --cl-config "max_table_rows: 10000" \\
    --cl-config "use_filename_as_sample_name: [{modules_conf}]" \\
    {extra_config_param}

    ls output/{report_filename}_data

    cp output/{report_filename}.html {j.html}
    cp output/{report_filename}_data/multiqc_data.json {j.json}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'

    if out_html_url:
        slack_message_cmd(j, text=f'*[{dataset_name}]* <{out_html_url}|MultiQC report>')

    j.command(wrap_command(cmd, setup_gcp=True))
    b.write_output(j.html, str(out_html_path))
    b.write_output(j.json, str(out_json_path))
    return j


def _write_sample_id_map(sample_map: dict[str, str], out_path: Path):
    """
    Configuring MultiQC to support bulk sample rename. `sample_map` is a dictionary
    of sample IDs. The map doesn't have to have records for all samples.
    Example:
    {
        'SID1': 'Patient1',
        'SID2': 'Patient2'
    }
    https://multiqc.info/docs/#bulk-sample-renaming-in-reports
    """
    with out_path.open('w') as fh:
        for sid, new_sid in sample_map.items():
            fh.write('\t'.join([sid, new_sid]) + '\n')
