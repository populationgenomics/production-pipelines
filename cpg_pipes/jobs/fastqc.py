"""
Jobs to run FastQC
"""
from os.path import basename
from typing import cast

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.providers.images import Images
from cpg_pipes.types import AlignmentInput, CramPath, FastqPath, FastqPairs
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.providers.refdata import RefData


def fastqc(
    b: hb.Batch,
    output_html_path: Path,
    output_zip_path: Path,
    alignment_input: AlignmentInput,
    refs: RefData,
    images: Images,
    subsample: bool = True,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Adds FastQC jobs. If the input is a set of fqs, runs FastQC on each fq file.
    """

    def _fastqc_one(jname_, input_path: CramPath | FastqPath):
        j = b.new_job(jname_, job_attrs)
        j.image(images.get('fastqc'))
        threads = STANDARD.set_resources(j, ncpu=16).get_nthreads()

        cmd = ''
        input_file = b.read_input(str(input_path))
        if isinstance(input_path, CramPath) and not input_path.is_bam:
            new_file = f'/io/batch/{basename(str(input_path))}'
            # FastQC doesn't support CRAMs, converting CRAM->BAM
            cmd += f"""\
            samtools view \\
            -T {refs.fasta_res_group(b).base} \\
            -@{threads - 1} \\
            -b {input_file} \\
            {"--subsample 0.01" if subsample else ''} \\
            -O {new_file}
            rm {input_file} 
            """
            input_file = new_file

        fname = basename(str(input_path))
        if not str(input_path).endswith('.gz'):
            cmd += f"""\
            mv {input_file} {input_file}.gz
            """
            input_file = f'{input_file}.gz'
            fname += '.gz'

        if subsample and not isinstance(input_path, CramPath):
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
        for lane_i, pair in enumerate(alignment_input):
            jname = f'FastQC R1'
            if len(alignment_input) > 1:
                jname += f' lane={lane_i}'
            j = _fastqc_one(jname, pair.r1)
            jobs.append(j)
    return jobs
