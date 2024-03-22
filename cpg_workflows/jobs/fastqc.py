"""
Jobs to run FastQC.
"""

from os.path import basename

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import command
from cpg_workflows.filetypes import (
    BamPath,
    FastqPath,
)
from cpg_workflows.resources import STANDARD


def fastqc(
    b: hb.Batch,
    output_html_path: Path,
    output_zip_path: Path,
    input_path: BamPath | FastqPath,
    subsample: bool = True,
    job_attrs: dict | None = None,
) -> Job:
    """
    Adds FastQC jobs. If the input is a set of fqs, runs FastQC on each fq file.
    """
    j = b.new_job('FASTQC', (job_attrs or {}) | {'tool': 'fastqc'})
    j.image(image_path('fastqc'))
    threads = STANDARD.set_resources(j, ncpu=16).get_nthreads()

    cmd = ''
    input_file: str | hb.ResourceFile = b.read_input(str(input_path))
    fname = basename(str(input_path))
    if not isinstance(input_path, BamPath):
        if not str(input_path).endswith('.gz'):
            cmd += f"""\
            mv {input_file} {input_file}.gz
            """
            input_file = f'{input_file}.gz'
            fname += '.gz'
        if subsample:
            n = 10_000_000
            new_file = f'$BATCH_TMPDIR/{fname}'
            cmd += f"""\
            gunzip -c {input_file} | head -n{n * 4} | gzip -c > {new_file}
            """
            input_file = new_file

    cmd += f"""\
    mkdir -p $BATCH_TMPDIR/outdir
    fastqc -t{threads} {input_file} \\
    --outdir $BATCH_TMPDIR/outdir
    ls $BATCH_TMPDIR/outdir
    ln $BATCH_TMPDIR/outdir/*_fastqc.html {j.out_html}
    ln $BATCH_TMPDIR/outdir/*_fastqc.zip {j.out_zip}
    unzip $BATCH_TMPDIR/outdir/*_fastqc.zip
    """
    j.command(command(cmd, monitor_space=True))
    b.write_output(j.out_html, str(output_html_path))
    b.write_output(j.out_zip, str(output_zip_path))
    return j
