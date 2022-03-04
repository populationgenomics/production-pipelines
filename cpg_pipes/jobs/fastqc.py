from os.path import join
from typing import Optional, List

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes.hb.command import wrap_command
from cpg_pipes.filetypes import AlignmentInput
from cpg_pipes.jobs import align


def fastqc(
    b: hb.Batch, 
    results_bucket: str,
    sample_name: str, 
    dataset_name: Optional[str], 
    alignment_input: AlignmentInput,
) -> List[Job]:
    """
    Adds FastQC jobs. If the input is a set of fqs, runs FastQC on each fq file.
    """
    def _fastqc_one(name, inp):
        j = b.new_job(name, dict(sample=sample_name, dataset=dataset_name))
        j.image('biocontainers/fastqc:v0.11.9_cv8')
        j.cpu(32)
        j.storage('150G')

        j.command(wrap_command(f"""\
        mkdir -p /io/batch/outdir
        fastqc -t64 {inp} --outdir /io/batch/outdir
        ls /io/batch/outdir
        ln /io/batch/outdir/*_fastqc.html {j.out_html}
        """, monitor_space=True))
    
        out_fname = name.replace(' ', '_').replace('/', '_').replace('=', '_') + '_fastqc.html'
        b.write_output(j.out_html, join(results_bucket, out_fname))
        return j

    jobs = []
    if alignment_input.bam_or_cram_path and alignment_input.bam_or_cram_path.endswith('.bam'):
        bam = alignment_input.as_cram_input_group(b)
        j = _fastqc_one('FastQC', bam.base)
        jobs.append(j)
        return jobs

    if alignment_input.bam_or_cram_path:
        assert alignment_input.bam_or_cram_path.endswith('.cram'), alignment_input
        extract_j = align.extract_fastq(
            b=b,
            cram=alignment_input.as_cram_input_group(b),
            sample_name=sample_name,
            dataset_name=dataset_name,
        )
        fqs1 = [extract_j.fq1]
        fqs2 = [extract_j.fq2]

    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        fqs1 = [b.read_input(fq) for fq in alignment_input.fqs1]
        fqs2 = [b.read_input(fq) for fq in alignment_input.fqs2]

    for lane_i in range(len(fqs1)):
        for r_i, fq in enumerate([fqs1[lane_i], fqs2[lane_i]]):
            name = f'FastQC R={r_i + 1}'
            if len(fqs1) > 1:
                name += f' lane={lane_i}'
            j = _fastqc_one(name, fq)
            jobs.append(j)
    return jobs
