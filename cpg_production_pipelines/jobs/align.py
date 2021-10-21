from enum import Enum
from textwrap import dedent, indent
from typing import Optional, List, Tuple
from os.path import splitext, basename, dirname, join
import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines import resources
from cpg_production_pipelines.jobs import wrap_command, new_job
from cpg_production_pipelines.jobs import picard
from cpg_production_pipelines.smdb import SMDB
from cpg_production_pipelines.utils import AlignmentInput

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
    output_path: str,
    sample_name: str,
    project_name: Optional[str] = None,
    aligner: Aligner = Aligner.BWA,
    markdup_tool: MarkDupTool = MarkDupTool.BIOBAMBAM,
    extra_label: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    overwrite: bool = False,
) -> Job:
    """
    - if the input is 1 fastq pair, submits one alignment job

    - if the input is a set of fastq pairs, submits multiple jobs per each pair
      then submits a separate merge job

    - if the input is a cram/bam
      - for bwa or bwa-mem2, stream bazam -> bwa
      - for dragmap, submit an extra job to extract a pair of fastqs from the cram/bam

    - if the markdup tool:
      - is biobambam2, stream the alignment or merging within the same job
      - is picard, submit a separate job with deduplication
    """
    ncpu = 32
    nthreads = ncpu * 2  # multithreading
    if alignment_input.fqs1 and len(alignment_input.fqs1) > 1:
        # running alignment in parallel, then merging
        assert len(alignment_input.fqs1) == len(alignment_input.fqs2), alignment_input
        align_jobs = []
        for i, (fq1, fq2) in enumerate(zip(alignment_input.fqs1, alignment_input.fqs2)):
            j, cmd = _align_one(
                b=b,
                alignment_input=AlignmentInput(fqs1=[fq1], fqs2=[fq2]),
                ncpu=ncpu,
                sample=sample_name,
                project=project_name,
                aligner=aligner,
                extra_label=f'{i+1}/{alignment_input.fqs1} {extra_label}',
            )
            cmd += sort_cmd(nthreads) + f' -o {j.sorted_bam}'
            j.command(wrap_command(cmd, monitor_space=True))
            align_jobs.append(j)
        first_j = align_jobs[0]

        merge_j = new_job('merge BAMs', sample_name, project_name)
        ncpu = 2
        nthreads = ncpu * 2  # multithreading
        merge_j.cpu(ncpu)
        merge_j.image(resources.SAMTOOLS_PICARD_IMAGE)
        # 300G for input BAMs, 300G for output BAM
        merge_j.storage('600G')

        align_cmd = f"""\
        samtools merge -@{nthreads-1} - {' '.join([j.sorted_bam for j in align_jobs])}
        """
        align_j = merge_j
        stdout_is_sorted = True

    else:
        align_j, align_cmd = _align_one(
            b=b,
            alignment_input=alignment_input,
            ncpu=ncpu,
            sample=sample_name,
            project=project_name,
            aligner=aligner,
            extra_label=extra_label,
        )
        first_j = align_j
        stdout_is_sorted = False

    last_j = finalise_alignment(
        b=b,
        align_cmd=align_cmd,
        stdout_is_sorted=stdout_is_sorted,
        j=align_j,
        sample_name=sample_name,
        project_name=project_name,
        nthreads=nthreads,
        output_path=output_path,
        overwrite=overwrite,
        markdup_tool=markdup_tool,
    )

    if depends_on:
        first_j.depends_on(*depends_on)
    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type='cram',
            output_path=output_path,
            sample_name=sample_name,
            project_name=project_name,
            first_j=first_j,
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def _align_one(
    b,
    alignment_input: AlignmentInput,
    ncpu: int,
    sample: str,
    project: Optional[str] = None,
    aligner: Aligner = Aligner.BWA,
    extra_label: Optional[str] = None,
) -> Tuple[Job, str]:
    """
    Creates a command that (re)aligns reads to hg38.

    It leaves sorting and duplicate marking to the user, thus returns a command in
    a raw string in addition to the Job object.
    """
    job_name = aligner.name
    if extra_label:
        job_name += f' {extra_label}'

    j = new_job(b, job_name, sample, project)
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.memory('standard')

    if aligner in [Aligner.BWAMEM2, Aligner.BWA]:
        if aligner == Aligner.BWAMEM2:
            tool = 'bwa-mem2'
            j.image(resources.BWAMEM2_IMAGE)
            index_exts = resources.BWAMEM2_INDEX_EXTS
        else:
            tool = 'bwa'
            j.image(resources.BWA_IMAGE)
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
            tool=tool,
        )
    else:
        j.image(resources.DRAGMAP_IMAGE)
        # Allow for 400G of extacted FQ, 150G of output CRAMs,
        # plus some tmp storage for sorting
        j.storage('700G')
        prep_inp_cmd = ''
        if alignment_input.bam_or_cram_path:
            extract_j = extract_fastq(
                b=b,
                cram=alignment_input.get_cram_input_group(b),
                sample_name=sample,
                project_name=project,
            )
            input_param = f'-1 {extract_j.fq1} -2 {extract_j.fq2}'

        else:
            # Allow for 300G input FQ, 150G output CRAM, plus some tmp storage
            if len(alignment_input.fqs1) == 1:
                input_param = (
                    f'-1 {b.read_input(alignment_input.fqs1[0])} '
                    f'-2 {b.read_input(alignment_input.fqs2[0])}'
                )
            else:
                files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
                files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
                prep_inp_cmd = f"""\
                cat {" ".join(files1)} >reads.R1.fq.gz
                cat {" ".join(files2)} >reads.R2.fq.gz
                # After merging lanes, we don't need input FQs anymore:
                rm {" ".join(files1 + files2)}  
                """
                input_param = f'-1 reads.R1.fq.gz -2 reads.R2.fq.gz'

        dragmap_index = b.read_input_group(
            **{
                k.replace('.', '_'): join(resources.DRAGMAP_INDEX_BUCKET, k)
                for k in resources.DRAGMAP_INDEX_FILES
            }
        )
        align_cmd = f"""\
        {indent(dedent(prep_inp_cmd), ' '*8)}
        dragen-os -r {dragmap_index} {input_param} \\
        --RGID {sample} --RGSM {sample}
        """

    return j, align_cmd


def _build_bwa_command(
    b,
    j,
    alignment_input,
    bwa_reference,
    nthreads,
    sample_name,
    tool,
) -> str:
    pull_inputs_cmd = ''
    if alignment_input.bam_or_cram_path:
        use_bazam = True
        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2

        if alignment_input.bam_or_cram_path.endswith('.cram'):
            storage = '400G'  # 150G input CRAM, 150G output CRAM, plus some tmp storage
        else:
            storage = '600G'  # input BAM is 1.5-2 times larger

        j.storage(storage)
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
        # Allow for 300G input FQ, 150G output CRAM, plus some tmp storage
        j.storage('600G')
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
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
    {tool} mem -K 100000000 {'-p' if use_bazam else ''} -t{nthreads} -Y \\
    -R '{rg_line}' {bwa_reference.base} {r1_param} {r2_param}
    """
    ).strip()


def extract_fastq(
    b: hb.Batch,
    cram: hb.ResourceGroup,
    sample_name: str,
    project_name: Optional[str] = None,
) -> Job:
    """
    Job that converts a bam or a cram to a pair of compressed fastq files
    """
    j = new_job(b, 'extract fastq', sample_name, project_name)
    ncpu = 32
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.image(resources.SAMTOOLS_PICARD_IMAGE)

    # 150G for input CRAM and 400G of extacted FASTQ
    j.storage('600G')

    reference = b.read_input_group(**resources.REF_D)
    cmd = f"""\
    samtools fastq -@{nthreads-1} {cram.base} -T {reference.base} \\
    -@{nthreads} -1 {j.fq1} -2 {j.fq2} -0 /dev/null -s /dev/null
    # After converting to FQ, we don't need input CRAMs anymore:
    rm {cram.base} {cram.index}
    """
    j.command(wrap_command(cmd, monitor_space=True))
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
    j.image(resources.DRAGMAP_IMAGE)
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
    | samtools sort -@{nthreads-1} -T /io/batch/samtools-sort-tmp -Obam
    """).strip()


def finalise_alignment(
    b: hb.Batch,
    align_cmd: str,
    stdout_is_sorted: bool,
    j: Job,
    sample_name: str,
    project_name: str,
    nthreads: int,
    output_path: Optional[str] = None,
    overwrite: bool = True,
    markdup_tool: Optional[MarkDupTool] = None,
) -> Tuple[str, Job]:

    reference = b.read_input_group(**resources.REF_D)

    md_j = None
    if markdup_tool == MarkDupTool.BIOBAMBAM:
        j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            }
        )
        align_cmd = dedent(f"""\
        {indent(dedent(align_cmd).strip(), ' '*8)} \\
        | bamsormadup inputformat=sam threads={nthreads} SO=coordinate \\
        M={j.duplicate_metrics} outputformat=sam \\
        tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp \\
        | samtools view -@{nthreads-1} -T {reference.base} -Ocram -o {j.output_cram.cram}       
        samtools index -@{nthreads-1} {j.output_cram.cram} {j.output_cram["cram.crai"]}
        """).strip()
        md_j = j
    else:
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
