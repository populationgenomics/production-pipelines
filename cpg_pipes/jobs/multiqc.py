#!/usr/bin/env python3

"""
Batch pipeline to run WGS QC.
"""

import logging
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.hb.batch import Batch
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.images import DRIVER_IMAGE

logger = logging.getLogger(__file__)


def multiqc(
    b: Batch,
    tmp_bucket: Path,
    paths: list[Path],
    out_html_path: Path,
    out_json_path: Path,
    out_html_url: str | None = None,
    ending_to_trim: set[str] | None = None,
    modules_to_trim_endings: set[str] | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run MultiQC for the files in `qc_paths`
    @param b: batch object
    @param tmp_bucket: bucket for tmp files
    @param paths: file bucket paths to pass into MultiQC 
    @param out_json_path: where to write MultiQC-generated JSON file
    @param out_html_path: where to write the HTML report
    @param out_html_url: URL corresponding to the HTML report
    @param ending_to_trim: trim these endings from input files to get sample names
    @param modules_to_trim_endings: list of modules for which trim the endings
    @param job_attrs: attributes to add to Hail Batch job
    @return: job object
    """
    j = b.new_job('Run MultiQC', job_attrs)
    j.image(DRIVER_IMAGE)
    STANDARD.set_resources(j, ncpu=16)

    file_list_path = tmp_bucket / 'multiqc-file-list.txt'
    with file_list_path.open('w') as f:
        f.writelines([f'{p}\n' for p in paths])
    file_list = b.read_input(str(file_list_path))

    endings_conf = ', '.join(list(ending_to_trim)) if ending_to_trim else ''
    modules_conf = ', '.join(
        list(modules_to_trim_endings)) if modules_to_trim_endings else ''

    cmd = f"""\
    pip install multiqc
    
    mkdir inputs
    cat {file_list} | gsutil -m cp -I inputs/

    multiqc inputs -o output -f \\
    --cl_config "extra_fn_clean_exts: [{endings_conf}]" \\
    --cl_config "use_filename_as_sample_name: [{modules_conf}]"

    cp output/multiqc_report.html {j.html}
    cp output/multiqc_data/multiqc_data.json {j.json}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'
    j.command(wrap_command(cmd, setup_gcp=True))

    b.write_output(j.html, str(out_html_path))
    b.write_output(j.json, str(out_json_path))
    return j
