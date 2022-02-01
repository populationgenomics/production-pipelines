"""
Create Hail Batch jobs for alignment.
"""
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
from cpg_pipes.hailbatch import AlignmentInput, PrevJob, wrap_command, \
    fasta_ref_resource, STANDARD

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


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


def samtools_stats(
    b,
    cram_path: str,
    sample_name: str,
    output_path: Optional[str] = None,
    project_name: Optional[str] = None,
    overwrite: bool = True,
    nthreads: Optional[int] = None,
) -> Job:
    """
    Run `samtools stats` for mapping QC
    """
    jname = 'samtools stats'
    j = b.new_job(jname, dict(sample=sample_name, project=project_name))
    if not output_path:
        output_path = cram_path + '.stats'
    if utils.can_reuse(cram_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.SAMTOOLS_PICARD_IMAGE)

    job_resource = STANDARD.set_resources(j, nthreads=nthreads)

    cram = b.read_input_group(**{
        'cram': cram_path,
        'cram.crai': cram_path + '.crai',
    })

    j.command(wrap_command(f"""\
    samtools stats -@{job_resource.get_nthreads() - 1} {cram.cram} > {j.output_stats}
    """))
    b.write_output(j.output_stats, output_path)
    
    return j


def align(
    b,
    alignment_input: AlignmentInput,
    sample_name: str,
    output_path: Optional[str] = None,
    project_name: Optional[str] = None,
    aligner: Aligner = Aligner.BWA,
    markdup_tool: MarkDupTool = MarkDupTool.BIOBAMBAM,
    extra_label: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    overwrite: bool = True,
    requested_nthreads: Optional[int] = None,
    number_of_shards_for_realignment: Optional[int] = None,
    prev_batch_jobs: Optional[Dict[Tuple[Optional[str], str], PrevJob]] = None,
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
    if output_path and utils.can_reuse(output_path, overwrite):
        job_name = aligner.name
        if extra_label:
            job_name += f' {extra_label}'
        job_name += ' [reuse]'
        return b.new_job(job_name, dict(sample=sample_name, project=project_name))

    if number_of_shards_for_realignment and number_of_shards_for_realignment > 1:
        if alignment_input.is_fastq():
            logger.warning(
                f'Cannot use number_of_shards_for_realignment for fastq inputs. '
                f'Sharding only works for CRAM/BAM inputs. '
                f'Sample: {project_name}/{sample_name}'
            )
            number_of_shards_for_realignment = None
    
    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_ncpu

    sharded_fq = alignment_input.is_fastq() and len(alignment_input.get_fqs1()) > 1
    sharded_bazam = (
        alignment_input.is_bam_or_cram() and 
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
            project=project_name,
            aligner=aligner,
        )
        first_j = align_j
        stdout_is_sorted = False
        output_fmt = 'sam'

    else:  # sharded alignment
        align_jobs = []
        sorted_bams = []

        if sharded_fq:
            # running alignment for each fastq pair in parallel
            fastq_pairs = zip(alignment_input.get_fqs1(), alignment_input.get_fqs2())
    
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
    
                # bwa-mem or dragmap command, but without sorting and deduplication:
                j, cmd = _align_one(
                    b=b,
                    job_name=jname,
                    alignment_input=AlignmentInput(fqs1=[fq1], fqs2=[fq2]),
                    requested_nthreads=requested_nthreads,
                    sample=sample_name,
                    project=project_name,
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
                    project=project_name,
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

        merge_j = b.new_job('Merge BAMs', dict(sample=sample_name, project=project_name))
        merge_j.image(resources.BIOINFO_IMAGE)
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
        project_name=project_name,
        requested_nthreads=requested_nthreads,
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
    job_name: str,
    alignment_input: AlignmentInput,
    requested_nthreads: int,
    sample: str,
    project: Optional[str] = None,
    aligner: Aligner = Aligner.BWA,
    number_of_shards_for_realignment: Optional[int] = None,
    shard_number_1based: Optional[int] = None,
 ) -> Tuple[Job, str]:
    """
    Creates a command that (re)aligns reads to hg38, and a Job object,
    but doesn't add the command to the Job object yet, so sorting and/or
    deduplication can be appended to the command.

    It leaves sorting and duplicate marking to the user, thus returns a command in
    a raw string in addition to the Job object.
    """
    j = b.new_job(job_name, dict(sample=sample, project=project, label=job_name))
    
    if number_of_shards_for_realignment is not None:
        assert number_of_shards_for_realignment > 1, number_of_shards_for_realignment
    
    if number_of_shards_for_realignment is None and requested_nthreads >= STANDARD.max_ncpu:
        # Running from FASTQ, or from CRAM without sharding, on a full instance.
        # We will need more storage, as default 350G (actual 344G) might not be enough, 
        # see example: https://batch.hail.populationgenomics.org.au/batches/7458/jobs/2
        job_resource = STANDARD.set_resources(j, attach_disk_storage_gb=550)
    else:
        job_resource = STANDARD.set_resources(j, nthreads=requested_nthreads)

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
            nthreads=job_resource.get_nthreads(),
            sample_name=sample,
            tool_name=tool_name,
            number_of_shards=number_of_shards_for_realignment,
            shard_number_1based=shard_number_1based,
        )
    else:
        j.image(resources.BIOINFO_IMAGE)
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
            input_param = f'-1 {extract_j.fq1} -2 {extract_j.fq2}'

        else:
            assert alignment_input.is_fastq()
            files1, files2 = alignment_input.as_fq_inputs(b)
            # Allow for 100G input FQ, 50G output CRAM, plus some tmp storage
            if len(alignment_input.get_fqs1()) == 1:
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
    number_of_shards: Optional[int] = None,
    shard_number_1based: Optional[int] = None,
) -> str:
    pull_inputs_cmd = ''
    if alignment_input.bam_or_cram_path:
        use_bazam = True
        assert not alignment_input.is_fastq()
        
        if alignment_input.bam_or_cram_path.startswith('gs://'):
            cram = alignment_input.as_cram_input_group(b)
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
        if number_of_shards and number_of_shards > 1:
            assert shard_number_1based is not None and shard_number_1based > 0, \
                (shard_number_1based, sample_name)
            shard_param = f' -s {shard_number_1based},{number_of_shards}'
        else:
            shard_param = ''
        bazam_cmd = (
            f'bazam -Xmx16g -Dsamjdk.reference_fasta={bwa_reference.base}'
            f' -n{min(nthreads, 6)} -bam {cram_localized_path}{shard_param} | '
        )
        r1_param = '-'
        r2_param = ''

    else:
        assert alignment_input.is_fastq()
        use_bazam = False
        bazam_cmd = ''
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
    {bazam_cmd} {tool_name} mem -K 100000000 {'-p' if use_bazam else ''} \\
    -t{nthreads - 1} -Y -R '{rg_line}' {bwa_reference.base} {r1_param} {r2_param}
    """).strip()


def extract_fastq(
    b,
    cram: hb.ResourceGroup,
    sample_name: str,
    project_name: Optional[str] = None,
    output_fq1: Optional[str] = None,
    output_fq2: Optional[str] = None,
) -> Job:
    """
    Job that converts a bam or a cram to an interleaved compressed fastq file
    """
    j = b.new_job('Extract fastq', dict(sample=sample_name, project=project_name))
    ncpu = 16
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.image(resources.BIOINFO_IMAGE)
    j.storage('700G')

    reference = fasta_ref_resource(b)
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
    reference = fasta_ref_resource(b)

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
    project_name: Optional[str],
    markdup_tool: MarkDupTool,
    output_path: Optional[str] = None,
    overwrite: bool = True,
    align_cmd_out_fmt: str = 'sam',
) -> Optional[Job]:
    """
    For MarkDupTool.BIOBAMBAM, adds bamsormadup command piped to the existing job.
    For MarkDupTool.PICARD, creates a new job, as Picard can't read from stdin.
    """

    reference = b.read_input_group(**resources.REF_D)
    
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
        | bamsormadup inputformat={align_cmd_out_fmt} threads={min(nthreads, 6)} SO=coordinate \\
        M={j.duplicate_metrics} outputformat=sam \\
        tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp \\
        | samtools view -@{min(nthreads, 6) - 1} -T {reference.base} -Ocram -o {j.output_cram.cram}       
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
