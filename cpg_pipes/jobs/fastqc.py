"""
Jobs to run FastQC
"""
from os.path import basename
from typing import cast

import hailtop.batch as hb
from cpg_utils.hail_batch import image_path
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.filetypes import (
    AlignmentInput,
    CramPath,
    FastqPath,
    FastqPairs,
    FastqPair,
)
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD


def fastqc(
    b: hb.Batch,
    output_html_path: Path,
    output_zip_path: Path,
    alignment_input: AlignmentInput,
    subsample: bool = True,
    job_attrs: dict | None = None,
    all_lanes: bool = False,
) -> list[Job]:
    """
    Adds FastQC jobs. If the input is a set of fqs, runs FastQC on each fq file.
    """
    if isinstance(alignment_input, CramPath) and not alignment_input.is_bam:
        raise NotImplementedError('FastQC does not support CRAM input')

    def _fastqc_one(jname_, input_path: CramPath | FastqPath):
        j = b.new_job(jname_, job_attrs)
        j.image(image_path('fastqc'))
        threads = STANDARD.set_resources(j, ncpu=16).get_nthreads()

        cmd = ''
        input_file = b.read_input(str(input_path))
        fname = basename(str(input_path))
        if isinstance(input_path, FastqPair):
            if not str(input_path).endswith('.gz'):
                cmd += f"""\
                mv {input_file} {input_file}.gz
                """
                input_file = f'{input_file}.gz'
                fname += '.gz'
            if subsample:
                n = 10_000_000
                new_file = f'/io/batch/{fname}'
                cmd += f"""\
                gunzip -c {input_file} | head -n{n * 4} | gzip -c > {new_file}
                """
                input_file = new_file

        cmd += f"""\
        mkdir -p /io/batch/outdir
        fastqc -t{threads} {input_file} \\
        --outdir /io/batch/outdir
        ls /io/batch/outdir
        ln /io/batch/outdir/*_fastqc.html {j.out_html}
        ln /io/batch/outdir/*_fastqc.zip {j.out_zip}
        unzip /io/batch/outdir/*_fastqc.zip
        """
        j.command(wrap_command(cmd, monitor_space=True))

        b.write_output(j.out_html, str(output_html_path))
        b.write_output(j.out_zip, str(output_zip_path))
        return j

    jobs = []
    if isinstance(alignment_input, CramPath):
        j = _fastqc_one('FastQC', cast(CramPath, alignment_input))
        jobs.append(j)
        return jobs
    else:
        assert isinstance(alignment_input, FastqPairs)
        if not all_lanes:
            alignment_input = FastqPairs([alignment_input[0]])
        for lane_i, pair in enumerate(alignment_input):
            jname = f'FastQC R1'
            if len(alignment_input) > 1:
                jname += f' lane={lane_i}'
            j = _fastqc_one(jname, pair.r1)
            jobs.append(j)
    return jobs
