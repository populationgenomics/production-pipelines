"""
Jobs to run FastQC
"""
from os.path import basename
from typing import cast

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.hail_batch import image_path
from cpg_utils.workflows.filetypes import (
    AlignmentInput,
    CramPath,
    FastqPath,
    FastqPairs,
    FastqPair,
)
from cpg_utils.hail_batch import command
from cpg_utils.workflows.resources import STANDARD


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
        j_ = b.new_job(jname_, job_attrs)
        j_.image(image_path('fastqc'))
        threads = STANDARD.set_resources(j_, ncpu=16).get_nthreads()

        cmd = ''
        input_file: str | hb.ResourceFile = b.read_input(str(input_path))
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
        ln $BATCH_TMPDIR/outdir/*_fastqc.html {j_.out_html}
        ln $BATCH_TMPDIR/outdir/*_fastqc.zip {j_.out_zip}
        unzip $BATCH_TMPDIR/outdir/*_fastqc.zip
        """
        j_.command(command(cmd, monitor_space=True))

        b.write_output(j_.out_html, str(output_html_path))
        b.write_output(j_.out_zip, str(output_zip_path))
        return j_

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
