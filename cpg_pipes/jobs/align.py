"""
Create Hail Batch jobs for alignment.
"""
from enum import Enum
from textwrap import dedent, indent
from typing import cast
from os.path import join
import logging

import hailtop.batch as hb
from cloudpathlib import CloudPath
from hailtop.batch.job import Job

from cpg_pipes import images, ref_data, buckets
from cpg_pipes.hb.prev_job import PrevJob
from cpg_pipes.jobs import picard
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.pipeline.analysis import AlignmentInput, AnalysisType, CramPath, \
    FastqPairs
from cpg_pipes.ref_data import REF_D

logger = logging.getLogger(__file__)


class Aligner(Enum):
    BWA = 'bwa'
    BWAMEM2 = 'bwa-mem2'
    DRAGMAP = 'dragmap'


class MarkDupTool(Enum):
    PICARD = 'picard'
    BIOBAMBAM = 'biobambam'
    NO_MARKDUP = 'no_markdup'


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
    output_path: str | CloudPath | None = None,
    qc_bucket: str | CloudPath | None = None,
    dataset_name: str | None = None,
    aligner: Aligner = Aligner.BWA,
    markdup_tool: MarkDupTool = MarkDupTool.BIOBAMBAM,
    extra_label: str | None = None,
    depends_on: list[Job] | None = None,
    smdb=None,
    overwrite: bool = True,
    requested_nthreads: int | None = None,
    number_of_shards_for_realignment: int | None = None,
    prev_batch_jobs: dict[tuple[str | None, str], PrevJob] | None = None,
) -> Job:
    """
    - if the input is 1 fastq pair, submits one alignment job.

    - if the input is a set of fastq pairs, submits multiple jobs per each pair,
      then submits a separate merge job.

    - if the input is a cram/bam
      - for bwa or bwa-mem2, stream bazam -> bwa.
        - if number_of_shards_for_realignment is > 1, use bazam to shard inputs
          and align in parallel, then merge result together.

      - for dragmap, submit an extra job to extract a pair of fastqs from the cram/bam,
        because dragmap can't read streamed files from bazam.

    - if the markdup tool:
      - is biobambam2, stream the alignment or merging within the same job.
      - is picard, submit a separate job with deduplication.
      
    - nthreads can be set for smaller test runs on toy instance, so the job 
      doesn't take entire 32-cpu/64-threaded instance.
    """
    if output_path and buckets.can_reuse(output_path, overwrite):
        job_name = aligner.name
        if extra_label:
            job_name += f' {extra_label}'
        job_name += ' [reuse]'
        return b.new_job(job_name, dict(sample=sample_name, dataset=dataset_name))

    if number_of_shards_for_realignment and number_of_shards_for_realignment > 1:
        if not isinstance(alignment_input, CramPath):
            logger.warning(
                f'Cannot use number_of_shards_for_realignment for fastq inputs. '
                f'Sharding only works for CRAM/BAM inputs. '
                f'Sample: {dataset_name}/{sample_name}'
            )
            number_of_shards_for_realignment = None
    
    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_threads()

    sharded_fq = (
        not isinstance(alignment_input, CramPath) and len(alignment_input) > 1
    )
    sharded_bazam = (
        isinstance(alignment_input, CramPath) and 
        number_of_shards_for_realignment and number_of_shards_for_realignment > 1
    )
    sharded = sharded_fq or sharded_bazam

    if not sharded:
        jname = f'{aligner.name}'
        if extra_label:
            jname += f' {extra_label}'

        align_j, align_cmd = _align_one(
            b=b,
            job_name=jname,
            alignment_input=alignment_input,
            requested_nthreads=requested_nthreads,
            sample=sample_name,
            dataset=dataset_name,
            aligner=aligner,
            storage_gb=500,
        )
        first_j = align_j
        stdout_is_sorted = False
        output_fmt = 'sam'

    else:  # sharded alignment
        align_jobs = []
        sorted_bams = []

        if sharded_fq:
            # running alignment for each fastq pair in parallel
            fastq_pairs = cast(FastqPairs, alignment_input)
            for i, pair in enumerate(fastq_pairs):
                jname = (
                    f'{aligner.name} {i+1}/{pair.r1}' + 
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
                    if buckets.can_reuse(existing_sorted_bam_path, overwrite):
                        logger.info(f'Reusing previous batch result: {existing_sorted_bam_path}')
                        jname += ' [reuse from previous batch]'
                        j = b.new_job(jname, dict(sample=sample_name, dataset=dataset_name))
                        align_jobs.append(j)
                        sorted_bams.append(b.read_input(existing_sorted_bam_path))
                        continue
    
                # bwa-mem or dragmap command, but without sorting and deduplication:
                j, cmd = _align_one(
                    b=b,
                    job_name=jname,
                    alignment_input=[pair],
                    requested_nthreads=requested_nthreads,
                    sample=sample_name,
                    dataset=dataset_name,
                    aligner=aligner,
                )
                cmd = cmd.strip()
                # Sorting with samtools, but not adding deduplication yet, because we
                # need to merge first.
                cmd += ' ' + sort_cmd(requested_nthreads) + f' -o {j.sorted_bam}'
                j.command(wrap_command(cmd, monitor_space=True))
                sorted_bams.append(j.sorted_bam)
                align_jobs.append(j)

        elif sharded_bazam:
            # running shared alignment for a CRAM, sharding with Bazam
            assert number_of_shards_for_realignment
            for shard_number_0based in range(number_of_shards_for_realignment):
                shard_number_1based = shard_number_0based + 1
                jname = (
                    f'{aligner.name} {shard_number_1based}/{number_of_shards_for_realignment}' + 
                    (f' {extra_label}' if extra_label else '')
                )
                j, cmd = _align_one(
                    b=b,
                    job_name=jname,
                    alignment_input=alignment_input,
                    sample=sample_name,
                    dataset=dataset_name,
                    aligner=aligner,
                    requested_nthreads=requested_nthreads,
                    number_of_shards_for_realignment=number_of_shards_for_realignment,
                    shard_number_1based=shard_number_1based,
                )
                cmd = cmd.strip()
                # Sorting with samtools, but not adding deduplication yet, because we
                # need to merge first.
                cmd += ' ' + sort_cmd(requested_nthreads) + f' -o {j.sorted_bam}'
                j.command(wrap_command(cmd, monitor_space=True))
                sorted_bams.append(j.sorted_bam)
                align_jobs.append(j)

        merge_j = b.new_job('Merge BAMs', dict(sample=sample_name, dataset=dataset_name))
        merge_j.image(images.BIOINFO_IMAGE)
        nthreads = STANDARD.set_resources(merge_j, nthreads=requested_nthreads).get_nthreads()

        align_cmd = f"""\
        samtools merge -@{nthreads - 1} - {' '.join(sorted_bams)}
        """.strip()
        output_fmt = 'bam'
        align_j = merge_j
        stdout_is_sorted = True
        first_j = merge_j
        
    md_j = finalise_alignment(
        b=b,
        align_cmd=align_cmd,
        stdout_is_sorted=stdout_is_sorted,
        j=align_j,
        sample_name=sample_name,
        dataset_name=dataset_name,
        requested_nthreads=requested_nthreads,
        markdup_tool=markdup_tool,
        output_path=output_path,
        qc_bucket=qc_bucket,
        overwrite=overwrite,
        align_cmd_out_fmt=output_fmt,
    )
    last_j = md_j if md_j is not None else align_j

    if depends_on:
        first_j.depends_on(*depends_on)
    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type=AnalysisType.CRAM,
            output_path=output_path,
            sample_names=[sample_name],
            dataset_name=dataset_name,
            first_j=first_j,
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def _align_one(
    b,
    job_name: str,
    alignment_input: AlignmentInput,
    requested_nthreads: int,
    sample: str,
    dataset: str | None = None,
    aligner: Aligner = Aligner.BWA,
    number_of_shards_for_realignment: int | None = None,
    shard_number_1based: int | None = None,
    storage_gb: float | None = None,
 ) -> tuple[Job, str]:
    """
    Creates a command that (re)aligns reads to hg38, and a Job object,
    but doesn't add the command to the Job object yet, so sorting and/or
    deduplication can be appended to the command.

    It leaves sorting and duplicate marking to the user, thus returns a command in
    a raw string in addition to the Job object.
    """
    j = b.new_job(job_name, dict(sample=sample, dataset=dataset, label=job_name))
    
    if number_of_shards_for_realignment is not None:
        assert number_of_shards_for_realignment > 1, number_of_shards_for_realignment

    job_resource = STANDARD.set_resources(
        j, nthreads=requested_nthreads, storage_gb=storage_gb
    )

    if aligner in [Aligner.BWAMEM2, Aligner.BWA]:
        if aligner == Aligner.BWAMEM2:
            tool_name = 'bwa-mem2'
            j.image(images.BWAMEM2_IMAGE)
            index_exts = ref_data.BWAMEM2_INDEX_EXTS
        else:
            tool_name = 'bwa'
            j.image(images.BIOINFO_IMAGE)
            index_exts = ref_data.BWA_INDEX_EXTS

        bwa_reference = b.read_input_group(
            **ref_data.REF_D,
            **{k: f'{ref_data.REF_FASTA}.{k}' for k in index_exts},
        )
        align_cmd = _build_bwa_command(
            b=b,
            alignment_input=alignment_input,
            bwa_reference=bwa_reference,
            nthreads=job_resource.get_nthreads(),
            sample_name=sample,
            tool_name=tool_name,
            number_of_shards=number_of_shards_for_realignment,
            shard_number_1based=shard_number_1based,
        )
    else:
        j.image(images.BIOINFO_IMAGE)
        dragmap_index = b.read_input_group(
            **{
                k.replace('.', '_'): join(ref_data.DRAGMAP_INDEX_BUCKET, k)
                for k in ref_data.DRAGMAP_INDEX_FILES
            }
        )
        prep_inp_cmd = ''
        if isinstance(alignment_input, CramPath):
            extract_j = extract_fastq(
                b=b,
                cram=alignment_input.resource_group(b),
                sample_name=sample,
                dataset_name=dataset,
            )
            input_param = f'-1 {extract_j.fq1} -2 {extract_j.fq2}'

        else:
            fastq_pairs = [p.as_resources(b) for p in cast(FastqPairs, alignment_input)]
            files1 = [pair[0] for pair in fastq_pairs]
            files2 = [pair[1] for pair in fastq_pairs]

            # Allow for 100G input FQ, 50G output CRAM, plus some tmp storage
            if len(fastq_pairs) == 1:
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
    alignment_input: AlignmentInput,
    bwa_reference: hb.ResourceGroup,
    nthreads: int,
    sample_name: str,
    tool_name: str,
    number_of_shards: int | None = None,
    shard_number_1based: int | None = None,
) -> str:
    pull_inputs_cmd = ''
    if isinstance(alignment_input, CramPath):
        use_bazam = True

        cram = cast(CramPath, alignment_input).resource_group(b)
        if number_of_shards and number_of_shards > 1:
            assert shard_number_1based is not None and shard_number_1based > 0, \
                (shard_number_1based, sample_name)
            shard_param = f' -s {shard_number_1based},{number_of_shards}'
        else:
            shard_param = ''
        bazam_cmd = (
            f'bazam -Xmx16g -Dsamjdk.reference_fasta={bwa_reference.base}'
            f' -n{min(nthreads, 6)} -bam {cram.cram}{shard_param} | '
        )
        r1_param = '-'
        r2_param = ''

    else:
        use_bazam = False
        bazam_cmd = ''
        fastq_pairs = [p.as_resources(b) for p in cast(FastqPairs, alignment_input)]
        files1 = [pair[0] for pair in fastq_pairs]
        files2 = [pair[1] for pair in fastq_pairs]
        if len(fastq_pairs) > 1:
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
    {bazam_cmd} {tool_name} mem -K 100000000 {'-p' if use_bazam else ''} \\
    -t{nthreads - 1} -Y -R '{rg_line}' {bwa_reference.base} {r1_param} {r2_param}
    """).strip()


def extract_fastq(
    b,
    cram: hb.ResourceGroup,
    sample_name: str,
    dataset_name: str | None = None,
    output_fq1: str | None = None,
    output_fq2: str | None = None,
) -> Job:
    """
    Job that converts a bam or a cram to an interleaved compressed fastq file
    """
    j = b.new_job('Extract fastq', dict(sample=sample_name, dataset=dataset_name))
    ncpu = 16
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.image(images.BIOINFO_IMAGE)
    j.storage('700G')

    reference = b.read_input_group(**REF_D)
    cmd = f"""\
    bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
    -n{nthreads} -bam {cram.base} -r1 {j.fq1} -r2 {j.fq2}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if output_fq1 or output_fq2:
        assert output_fq1 and output_fq2, (output_fq1, output_fq2)
        b.write_output(j.fq1, output_fq1)
        b.write_output(j.fq2, output_fq2)
    return j


def create_dragmap_index(b: hb.Batch) -> Job:
    """
    Creates the index for DRAGMAP
    """
    reference = b.read_input_group(**REF_D)

    j = b.new_job('Index DRAGMAP')
    j.image(images.BIOINFO_IMAGE)
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
    for f in ref_data.DRAGMAP_INDEX_FILES:
        cmd += f'ln $DIR/{f} {getattr(j, f.replace(".", "_"))}\n'
    cmd += 'df -h; pwd; ls | grep -v proc | xargs du -sh'
    j.command(cmd)
    for f in ref_data.DRAGMAP_INDEX_FILES:
        b.write_output(
            getattr(j, f.replace('.', '_')), join(ref_data.DRAGMAP_INDEX_BUCKET, f)
        )
    return j


def sort_cmd(requested_nthreads: int):
    nthreads = STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()
    return dedent(f"""\
    | samtools sort -@{min(nthreads, 6) - 1} -T /io/batch/samtools-sort-tmp -Obam
    """).strip()


def finalise_alignment(
    b: hb.Batch,
    align_cmd: str,
    stdout_is_sorted: bool,
    j: Job,
    requested_nthreads: int,
    sample_name: str,
    dataset_name: str | None,
    markdup_tool: MarkDupTool,
    output_path: CloudPath | str | None = None,
    qc_bucket: CloudPath | str | None = None,
    overwrite: bool = True,
    align_cmd_out_fmt: str = 'sam',
) -> Job | None:
    """
    For MarkDupTool.BIOBAMBAM, adds bamsormadup command piped to the existing job.
    For MarkDupTool.PICARD, creates a new job, as Picard can't read from stdin.
    """

    reference = b.read_input_group(**ref_data.REF_D)
    
    nthreads = STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()

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
        | bamsormadup inputformat={align_cmd_out_fmt} threads={min(nthreads, 6)} \\
        SO=coordinate M={j.duplicate_metrics} outputformat=sam \\
        tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp \\
        | samtools view -@{min(nthreads, 6) - 1} -T {reference.base} \\
        -Ocram -o {j.output_cram.cram}       
        
        samtools index -@{nthreads - 1} {j.output_cram.cram} \\
        {j.output_cram["cram.crai"]}
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
            dataset_name=dataset_name,
            overwrite=overwrite,
        )

    if output_path:
        output_path = CloudPath(output_path)

        if qc_bucket:
            qc_bucket = CloudPath(qc_bucket)
        else:
            qc_bucket = output_path.parent
    
        if md_j is not None:
            b.write_output(md_j.output_cram, str(output_path.with_suffix('')))
            b.write_output(md_j.duplicate_metrics, str(
                qc_bucket /
                'duplicate-metrics' /
                f'{sample_name}-duplicate-metrics.csv'
            ))
        else:
            b.write_output(j.sorted_bam, str(output_path.with_suffix('')))

    return md_j
