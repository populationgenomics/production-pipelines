"""
Create Hail Batch jobs for alignment.
"""
import math
from enum import Enum
from textwrap import dedent, indent
from typing import Optional, List, Tuple, Dict
from os.path import splitext, basename, dirname, join
import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import resources, utils
from cpg_pipes.jobs import picard
from cpg_pipes.smdb import SMDB
from cpg_pipes.hailbatch import AlignmentInput, PrevJob, wrap_command

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Aligner(Enum):
    BWA = 1
    BWAMEM2 = 2
    DRAGMAP = 3


class MarkDupTool(Enum):
    PICARD = 1
    BIOBAMBAM = 2
    NO_MARKDUP = 3


def dragmap(*args, **kwargs):
    """
    Runs DRAGMAP to (re)align reads to hg38, and outputs sorted CRAM with
    marked duplicates
    """
    kwargs['aligner'] = Aligner.DRAGMAP
    return align(*args, **kwargs)


def bwa(*args, **kwargs):
    """
    Runs BWA to (re)align reads to hg38, and outputs sorted CRAM with
    marked duplicates
    """
    kwargs['aligner'] = Aligner.BWA
    return align(*args, **kwargs)


def bwamem2(*args, **kwargs):
    """
    Runs BWA-MEM2 to (re)align reads to hg38, and outputs sorted CRAM with
    marked duplicates
    """
    kwargs['aligner'] = Aligner.BWAMEM2
    return align(*args, **kwargs)


def align(
    b,
    alignment_input: AlignmentInput,
    sample_name: str,
    output_path: Optional[str] = None,
    project_name: Optional[str] = None,
    aligner: Aligner = Aligner.DRAGMAP,
    markdup_tool: MarkDupTool = MarkDupTool.PICARD,
    extra_label: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    overwrite: bool = True,
    check_existence: bool = True,
    prev_batch_jobs: Optional[Dict[Tuple[Optional[str], str], PrevJob]] = None,
    nthreads: Optional[int] = None,
) -> Job:
    """
    - if the input is 1 fastq pair, submits one alignment job

    - if the input is a set of fastq pairs, submits multiple jobs per each pair
      then submits a separate merge job

    - if the input is a cram/bam
      - for bwa or bwa-mem2, stream bazam -> bwa
      - for dragmap, submit an extra job to extract a pair of fastqs from the cram/bam,
        because dragmap can't read streamed files from bazam

    - if the markdup tool:
      - is biobambam2, stream the alignment or merging within the same job
      - is picard, submit a separate job with deduplication
      
    - fraction_of_64thread_instance can be set for smaller test runs on toy instance,
      so the job doesn't acquire a whole 32-cpu/64-threaded instance
    """
    if output_path and check_existence and utils.can_reuse(output_path, overwrite):
        job_name = aligner.name
        if extra_label:
            job_name += f' {extra_label}'
        job_name += ' [reuse]'
        return b.new_job(job_name, dict(sample=sample_name, project=project_name))

    total_ncpu = 16
    storage_gb = 375
    if nthreads is not None:
        if nthreads < 1:
            raise ValueError('align: nthreads must be >= 1')
        # round to the nearest power of 2 (15 -> 16, 16 -> 16, 17 -> 32)
        nthreads = int(pow(2, math.ceil(math.log2(nthreads))))
        ncpu = nthreads * 2  # multithreading
        # When requesting a 16- or 32-cpu job, the job will occupy entire instance
        # with a fixed 375GB SSD. When we want a n'th fraction of cpus, Batch will 
        # squeeze N jobs like that on an instance, thus the storage will be split 
        # N-way accordingly.
        fraction = ncpu / total_ncpu
        storage_gb = int(math.ceil(storage_gb * fraction))
    else:
        ncpu = total_ncpu
        nthreads = ncpu * 2

    if alignment_input.fqs1 is not None and len(alignment_input.fqs1) > 1:
        # running alignment in parallel, then merging
        assert alignment_input.fqs2 is not None
        assert len(alignment_input.fqs1) == len(alignment_input.fqs2), alignment_input
        fastq_pairs = zip(alignment_input.fqs1, alignment_input.fqs2)

        align_jobs = []
        sorted_bams = []
        for i, (fq1, fq2) in enumerate(fastq_pairs):
            jname = (
                f'{aligner.name} {i+1}/{fq1}' + 
                (f' {extra_label}' if extra_label else '')
            )
            key = sample_name, jname 
            if prev_batch_jobs and key in prev_batch_jobs:
                prevj = prev_batch_jobs[key]
                logger.info(
                    f'Job was run in the previous batch {prevj.batch_number}: {key}'
                )
                existing_sorted_bam_path = (
                    f'{prevj.hail_bucket}/batch/{prevj.batchid}/{prevj.job_number}/sorted_bam'
                )
                if utils.can_reuse(existing_sorted_bam_path, overwrite):
                    logger.info(f'Reusing previous batch result: {existing_sorted_bam_path}')
                    jname += ' [reuse from previous batch]'
                    j = b.new_job(jname, dict(sample=sample_name, project=project_name))
                    align_jobs.append(j)
                    sorted_bams.append(b.read_input(existing_sorted_bam_path))
                    continue
            j, cmd = _align_one(
                b=b,
                alignment_input=AlignmentInput(fqs1=[fq1], fqs2=[fq2]),
                job_name=jname,
                ncpu=ncpu,
                sample=sample_name,
                project=project_name,
                aligner=aligner,
            )
            j.storage(f'{storage_gb}G')
            cmd = cmd.strip()
            cmd += ' ' + sort_cmd(nthreads) + f' -o {j.sorted_bam}'
            j.command(wrap_command(cmd, monitor_space=True))
            sorted_bams.append(j.sorted_bam)
            align_jobs.append(j)

        merge_j = b.new_job('Merge BAMs', dict(sample=sample_name, project=project_name))
        merge_j.cpu(ncpu)
        merge_j.image(resources.BIOINFO_IMAGE)
        merge_j.storage(f'{storage_gb}G')

        align_cmd = f"""\
        samtools merge -@{nthreads - 1} - {' '.join(sorted_bams)}
        """.strip()
        output_fmt = 'bam'
        align_j = merge_j
        stdout_is_sorted = True
        # first_j = align_jobs[0]
        first_j = merge_j

    else:
        jname = f'{aligner.name}'
        if extra_label:
            jname += f' {extra_label}'
        align_j, align_cmd = _align_one(
            b=b,
            job_name=jname,
            alignment_input=alignment_input,
            ncpu=ncpu,
            sample=sample_name,
            project=project_name,
            aligner=aligner,
        )
        align_j.storage(f'{storage_gb}G')

        first_j = align_j
        stdout_is_sorted = False
        output_fmt = 'sam'

    md_j = finalise_alignment(
        b=b,
        align_cmd=align_cmd,
        stdout_is_sorted=stdout_is_sorted,
        j=align_j,
        sample_name=sample_name,
        project_name=project_name,
        nthreads=nthreads,
        markdup_tool=markdup_tool,
        output_path=output_path,
        overwrite=overwrite,
        align_cmd_out_fmt=output_fmt,
    )
    last_j = md_j if md_j is not None else align_j

    if depends_on:
        first_j.depends_on(*depends_on)
    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type='cram',
            output_path=output_path,
            sample_names=[sample_name],
            project_name=project_name,
            first_j=first_j,
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def _align_one(
    b,
    alignment_input: AlignmentInput,
    job_name: str,
    ncpu: int,
    sample: str,
    project: Optional[str] = None,
    aligner: Aligner = Aligner.BWA,
) -> Tuple[Job, str]:
    """
    Creates a command that (re)aligns reads to hg38.

    It leaves sorting and duplicate marking to the user, thus returns a command in
    a raw string in addition to the Job object.
    """
    j = b.new_job(job_name, dict(sample=sample, project=project, label=job_name))
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.memory('standard')

    if aligner in [Aligner.BWAMEM2, Aligner.BWA]:
        if aligner == Aligner.BWAMEM2:
            tool_name = 'bwa-mem2'
            j.image(resources.BWAMEM2_IMAGE)
            index_exts = resources.BWAMEM2_INDEX_EXTS
        else:
            tool_name = 'bwa'
            j.image(resources.BIOINFO_IMAGE)
            index_exts = resources.BWA_INDEX_EXTS

        bwa_reference = b.read_input_group(
            **resources.REF_D,
            **{k: f'{resources.REF_FASTA}.{k}' for k in index_exts},
        )
        align_cmd = _build_bwa_command(
            b=b,
            j=j,
            alignment_input=alignment_input,
            bwa_reference=bwa_reference,
            nthreads=nthreads,
            sample_name=sample,
            tool_name=tool_name,
        )
    else:
        j.image(resources.BIOINFO_IMAGE)
        j.storage(f'300G')
        dragmap_index = b.read_input_group(
            **{
                k.replace('.', '_'): join(resources.DRAGMAP_INDEX_BUCKET, k)
                for k in resources.DRAGMAP_INDEX_FILES
            }
        )
        prep_inp_cmd = ''
        if alignment_input.bam_or_cram_path:
            extract_j = extract_fastq(
                b=b,
                cram=alignment_input.as_cram_input_group(b),
                sample_name=sample,
                project_name=project,
            )
            input_param = f'-1 {extract_j.output_interleaved_fq} --interleaved=1'

        else:
            files1, files2 = alignment_input.as_fq_inputs(b)
            # Allow for 100G input FQ, 50G output CRAM, plus some tmp storage
            if alignment_input.fqs1 is not None and len(alignment_input.fqs1) == 1:
                input_param = f'-1 {files1[0]} -2 {files2[0]}'
            else:
                prep_inp_cmd = f"""\
                cat {" ".join(files1)} >reads.R1.fq.gz
                cat {" ".join(files2)} >reads.R2.fq.gz
                # After merging lanes, we don't need input FQs anymore:
                rm {" ".join(files1 + files2)}  
                """
                input_param = f'-1 reads.R1.fq.gz -2 reads.R2.fq.gz'

        align_cmd = f"""\
        {indent(dedent(prep_inp_cmd), ' '*8)}
        dragen-os -r {dragmap_index} {input_param} \\
        --RGID {sample} --RGSM {sample}
        """

    return j, align_cmd


def _build_bwa_command(
    b,
    j: Job,
    alignment_input: AlignmentInput,
    bwa_reference: hb.ResourceGroup,
    nthreads: int,
    sample_name: str,
    tool_name: str,
) -> str:
    pull_inputs_cmd = ''
    if alignment_input.bam_or_cram_path:
        use_bazam = True
        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        
        if alignment_input.bam_or_cram_path.startswith('gs://'):
            cram = b.read_input_group(
                base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
            )
            cram_localized_path = cram.base
        else:
            # Can't use on Batch localization mechanism with `b.read_input_group`,
            # but have to manually localize with `wget`
            cram_name = basename(alignment_input.bam_or_cram_path)
            work_dir = dirname(j.output_cram.cram)
            cram_localized_path = join(work_dir, cram_name)
            index_ext = '.crai' if cram_name.endswith('.cram') else '.bai'
            crai_localized_path = join(work_dir, cram_name + index_ext)
            pull_inputs_cmd = (
                f'wget {alignment_input.bam_or_cram_path} -O {cram_localized_path}\n'
                f'wget {alignment_input.index_path} -O {crai_localized_path}'
            )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={bwa_reference.base}'
            f' -n{nthreads} -bam {cram_localized_path})'
        )
        r2_param = '-'

    else:
        assert alignment_input.fqs1 and alignment_input.fqs2, alignment_input
        assert len(alignment_input.fqs1) == len(alignment_input.fqs1), alignment_input
        use_bazam = False
        files1, files2 = alignment_input.as_fq_inputs(b)
        if len(files1) > 1:
            r1_param = f'<(cat {" ".join(files1)})'
            r2_param = f'<(cat {" ".join(files2)})'
        else:
            r1_param = files1[0]
            r2_param = files2[0]

    rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
    # BWA command options:
    # -K   process INT input bases in each batch regardless of nThreads (for reproducibility)
    # -p   smart pairing (ignoring in2.fq)
    # -t16 threads
    # -Y   use soft clipping for supplementary alignments
    # -R   read group header line such as '@RG\tID:foo\tSM:bar'
    return dedent(f"""\
    {pull_inputs_cmd}
    {tool_name} mem -K 100000000 {'-p' if use_bazam else ''} -t{nthreads} -Y \\
    -R '{rg_line}' {bwa_reference.base} {r1_param} {r2_param}
    """).strip()


def extract_fastq(
    b,
    cram: hb.ResourceGroup,
    sample_name: str,
    project_name: Optional[str] = None,
    output_interleaved_fq: Optional[str] = None,
) -> Job:
    """
    Job that converts a bam or a cram to an interleaved compressed fastq file
    """
    j = b.new_job('Extract fastq', dict(sample=sample_name, project=project_name))
    ncpu = 16
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.image(resources.BIOINFO_IMAGE)
    j.storage('375G')

    reference = b.read_input_group(**resources.REF_D)
    cmd = f"""\
    bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
    -n{nthreads} -bam {cram.base} | gzip -c > {j.output_interleaved_fq}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if output_interleaved_fq:
        b.write_output(j.output_interleaved_fq, output_interleaved_fq)
    return j


def create_dragmap_index(b: hb.Batch) -> Job:
    """
    Creates the index for DRAGMAP
    """
    reference = b.read_input_group(
        base=resources.REF_FASTA,
        fai=resources.REF_FASTA + '.fai',
        dict=resources.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
    )
    j = b.new_job('Index DRAGMAP')
    j.image(resources.BIOINFO_IMAGE)
    j.memory('standard')
    j.cpu(32)
    j.storage('40G')
    cmd = f"""\
    DIR=$(dirname {j.hash_table_cfg})

    dragen-os \\
    --build-hash-table true \\
    --ht-reference {reference.base} \\
    --ht-num-threads 32 \\
    --output-directory $DIR
    """
    j.command(wrap_command(cmd))
    for f in resources.DRAGMAP_INDEX_FILES:
        cmd += f'ln $DIR/{f} {getattr(j, f.replace(".", "_"))}\n'
    cmd += 'df -h; pwd; ls | grep -v proc | xargs du -sh'
    j.command(cmd)
    for f in resources.DRAGMAP_INDEX_FILES:
        b.write_output(
            getattr(j, f.replace('.', '_')), join(resources.DRAGMAP_INDEX_BUCKET, f)
        )
    return j


def sort_cmd(nthreads):
    return dedent(f"""\
    | samtools sort -@{nthreads - 1} -T /io/batch/samtools-sort-tmp -Obam
    """).strip()


def finalise_alignment(
    b: hb.Batch,
    align_cmd: str,
    stdout_is_sorted: bool,
    j: Job,
    sample_name: str,
    project_name: Optional[str],
    nthreads: int,
    markdup_tool: MarkDupTool,
    output_path: Optional[str] = None,
    overwrite: bool = True,
    align_cmd_out_fmt: str = 'sam',
) -> Optional[Job]:

    reference = b.read_input_group(**resources.REF_D)

    md_j = None
    if markdup_tool == MarkDupTool.BIOBAMBAM:
        j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            }
        )
        align_cmd = f"""\
        {align_cmd.strip()} \\
        | bamsormadup inputformat={align_cmd_out_fmt} threads={nthreads} SO=coordinate \\
        M={j.duplicate_metrics} outputformat=sam \\
        tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp \\
        | samtools view -@{nthreads - 1} -T {reference.base} -Ocram -o {j.output_cram.cram}       
        samtools index -@{nthreads - 1} {j.output_cram.cram} {j.output_cram["cram.crai"]}
        """.strip()
        md_j = j
    else:
        align_cmd = align_cmd.strip()
        if not stdout_is_sorted:
            align_cmd += f' {sort_cmd(nthreads)}'
        align_cmd += f' > {j.sorted_bam}'

    j.command(wrap_command(align_cmd, monitor_space=True))

    if markdup_tool == MarkDupTool.PICARD:
        md_j = picard.markdup(
            b,
            j.sorted_bam,
            sample_name=sample_name,
            project_name=project_name,
            overwrite=overwrite,
        )

    if output_path:
        if md_j is not None:
            b.write_output(md_j.output_cram, splitext(output_path)[0])
            b.write_output(
                md_j.duplicate_metrics,
                join(
                    dirname(output_path),
                    'duplicate-metrics',
                    f'{sample_name}-duplicate-metrics.csv',
                ),
            )
        else:
            b.write_output(j.sorted_bam, splitext(output_path)[0])

    return md_j
