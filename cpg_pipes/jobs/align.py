"""
Create Hail Batch jobs for alignment.
"""

from enum import Enum
from textwrap import dedent
from typing import cast
import logging

import hailtop.batch as hb
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.targets import Sample
from cpg_pipes.filetypes import AlignmentInput, FastqPairs, CramPath
from cpg_pipes.jobs import picard
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.utils import can_reuse, exists

logger = logging.getLogger(__file__)


BWA_INDEX_EXTS = ['sa', 'amb', 'bwt', 'ann', 'pac', 'alt']
BWAMEM2_INDEX_EXTS = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac', 'alt']
DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']

DEFAULT_REALIGNMENT_SHARD_NUM = 10


class Aligner(Enum):
    """
    Tool that performs the alignment. Value must be the name of the executable.
    """

    BWA = 'bwa'
    BWAMEM2 = 'bwa-mem2'
    DRAGMAP = 'dragmap'


class MarkDupTool(Enum):
    """
    Tool that performs de-duplication.
    """

    PICARD = 'picard'
    BIOBAMBAM = 'biobambam'
    NO_MARKDUP = 'no_markdup'


def _get_cram_reference_from_version(cram_version) -> str:
    """
    Get the reference used for the specific cram_version,
    so that bazam is able to correctly decompress the reads
    """
    cram_version_map = get_config()['workflow'].get('cram_version_reference', {})
    if cram_version in cram_version_map:
        return cram_version_map[cram_version]
    raise ValueError(
        f'Unrecognised cram_version: "{cram_version}", expected one of: {", ".join(cram_version_map.keys())}'
    )


class MissingAlignmentInputException(Exception):
    """Raise if alignment input is missing"""

    pass


def _get_alignment_input(sample: Sample) -> AlignmentInput:
    sequencing_type = get_config()['workflow']['sequencing_type']
    alignment_input = sample.alignment_input_by_seq_type.get(sequencing_type)
    if realign_cram_ver := get_config()['workflow'].get('realign_from_cram_version'):
        if (
            path := (
                sample.dataset.prefix()
                / 'cram'
                / realign_cram_ver
                / f'{sample.id}.cram'
            )
        ).exists():
            logger.info(f'Realigning from {realign_cram_ver} CRAM {path}')
            alignment_input = CramPath(
                path,
                reference_assembly=_get_cram_reference_from_version(realign_cram_ver),
            )

    if alignment_input is None or not alignment_input.exists():
        raise MissingAlignmentInputException(
            f'No alignment inputs found for sample {sample}'
            + (f': {alignment_input}' if alignment_input else '')
        )
    return alignment_input


def align(
    b,
    sample: Sample,
    job_attrs: dict | None = None,
    output_path: CramPath | None = None,
    out_markdup_metrics_path: Path | None = None,
    aligner: Aligner = Aligner.DRAGMAP,
    markdup_tool: MarkDupTool = MarkDupTool.PICARD,
    extra_label: str | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
    realignment_shards_num: int = DEFAULT_REALIGNMENT_SHARD_NUM,
) -> list[Job]:
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
    if output_path and can_reuse(output_path.path, overwrite):
        return []

    alignment_input = _get_alignment_input(sample)

    base_jname = 'Align'
    if extra_label:
        base_jname += f' {extra_label}'

    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_threads()

    sharded_fq = isinstance(alignment_input, FastqPairs) and len(alignment_input) > 1
    sharded_bazam = isinstance(alignment_input, CramPath) and realignment_shards_num > 1
    sharded = sharded_fq or sharded_bazam

    jobs = []
    sharded_align_jobs = []
    sorted_bams = []

    if not sharded:  # Just running one alignment job
        align_j, align_cmd = _align_one(
            b=b,
            job_name=base_jname,
            alignment_input=alignment_input,
            requested_nthreads=requested_nthreads,
            sample_name=sample.id,
            job_attrs=job_attrs,
            aligner=aligner,
            should_sort=False,
        )
        stdout_is_sorted = False
        output_fmt = 'sam'
        jobs.append(align_j)
        merge_or_align_j = align_j

    else:  # Aligning in parallel and merging afterwards
        if sharded_fq:  # Aligning each lane separately, merging after
            # running alignment for each fastq pair in parallel
            fastq_pairs = cast(FastqPairs, alignment_input)
            for pair in fastq_pairs:
                # bwa-mem or dragmap command, but without sorting and deduplication:
                j, cmd = _align_one(
                    b=b,
                    job_name=base_jname,
                    alignment_input=FastqPairs([pair]),
                    requested_nthreads=requested_nthreads,
                    sample_name=sample.id,
                    job_attrs=job_attrs,
                    aligner=aligner,
                    should_sort=True,
                )
                j.command(wrap_command(cmd, monitor_space=True))
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        elif sharded_bazam:  # Using BAZAM to shard CRAM
            assert realignment_shards_num, realignment_shards_num
            for shard_number in range(realignment_shards_num):
                j, cmd = _align_one(
                    b=b,
                    job_name=base_jname,
                    alignment_input=alignment_input,
                    sample_name=sample.id,
                    job_attrs=job_attrs,
                    aligner=aligner,
                    requested_nthreads=requested_nthreads,
                    number_of_shards_for_realignment=realignment_shards_num,
                    shard_number=shard_number,
                    should_sort=True,
                )
                # Sorting with samtools, but not adding deduplication yet, because we
                # need to merge first.
                j.command(wrap_command(cmd, monitor_space=True))
                sorted_bams.append(str(j.sorted_bam))
                sharded_align_jobs.append(j)

        merge_j = b.new_job(
            'Merge BAMs', (job_attrs or {}) | dict(tool='samtools_merge')
        )
        merge_j.image(image_path('bwa'))

        nthreads = STANDARD.set_resources(
            merge_j,
            nthreads=requested_nthreads,
            # for FASTQ or BAM inputs, requesting more disk (400G). Example when
            # default is not enough: https://batch.hail.populationgenomics.org.au/batches/73892/jobs/56
            storage_gb=storage_for_align_job(
                alignment_input=alignment_input,
            ),
        ).get_nthreads()

        align_cmd = f"""\
        samtools merge -@{nthreads - 1} - {' '.join(sorted_bams)}
        """.strip()
        output_fmt = 'bam'
        jobs.extend(sharded_align_jobs)
        jobs.append(merge_j)
        merge_or_align_j = merge_j
        stdout_is_sorted = True

    md_j = finalise_alignment(
        b=b,
        align_cmd=align_cmd,
        stdout_is_sorted=stdout_is_sorted,
        j=merge_or_align_j,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
        markdup_tool=markdup_tool,
        output_path=output_path,
        out_markdup_metrics_path=out_markdup_metrics_path,
        align_cmd_out_fmt=output_fmt,
        overwrite=overwrite,
    )
    if md_j != merge_or_align_j:
        jobs.append(md_j)

    return jobs


def storage_for_align_job(alignment_input: AlignmentInput) -> int | None:
    """
    Get storage for an alignment job, gb
    """
    storage_gb = None  # avoid attaching extra disk by default

    try:
        storage_gb = get_config()['workflow']['resources']['Align']['storage_gb']
    except KeyError:
        pass
    else:
        assert isinstance(storage_gb, int), storage_gb
        return storage_gb

    sequencing_type = get_config()['workflow']['sequencing_type']
    if sequencing_type == 'genome':
        if (
            isinstance(alignment_input, FastqPairs)
            or isinstance(alignment_input, CramPath)
            and alignment_input.is_bam
        ):
            storage_gb = 400  # for WGS FASTQ or BAM inputs, need more disk
    return storage_gb


def _align_one(
    b,
    job_name: str,
    alignment_input: AlignmentInput,
    requested_nthreads: int,
    sample_name: str,
    job_attrs: dict | None = None,
    aligner: Aligner = Aligner.BWA,
    number_of_shards_for_realignment: int | None = None,
    shard_number: int | None = None,
    should_sort: bool = False,
) -> tuple[Job, str]:
    """
    Creates a job that (re)aligns reads to hg38. Returns the job object and a command
    separately, and doesn't add the command to the Job object, so stream-sorting
    and/or deduplication can be appended to the command later.
    """

    if number_of_shards_for_realignment is not None:
        assert number_of_shards_for_realignment > 1, number_of_shards_for_realignment

    job_attrs = (job_attrs or {}) | dict(label=job_name, tool=aligner.name)
    if shard_number is not None and number_of_shards_for_realignment is not None:
        job_name = (
            f'{job_name} ' f'{shard_number + 1}/{number_of_shards_for_realignment} '
        )
    job_name = f'{job_name} {alignment_input.path_glob()}'
    j = b.new_job(job_name, job_attrs)

    nthreads = STANDARD.set_resources(
        j,
        nthreads=requested_nthreads,
        storage_gb=storage_for_align_job(alignment_input=alignment_input),
    ).get_nthreads()

    # 2022-07-22 mfranklin:
    #   Replace process substitution with named-pipes (FIFO)
    #   This is named-pipe name -> command to populate it
    fifo_commands: dict[str, str] = {}

    index_cmd = ''
    if isinstance(alignment_input, CramPath):
        use_bazam = True
        if number_of_shards_for_realignment and number_of_shards_for_realignment > 1:
            assert shard_number is not None and shard_number >= 0, (
                shard_number,
                sample_name,
            )
            shard_param = f' -s {shard_number + 1},{number_of_shards_for_realignment}'
        else:
            shard_param = ''

        reference_inp = None
        if not alignment_input.is_bam:  # if is CRAM
            assert (
                alignment_input.reference_assembly
            ), f'The reference input for the alignment input "{alignment_input.path}" was not set'
            reference_inp = b.read_input_group(
                base=str(alignment_input.reference_assembly),
                fai=str(alignment_input.reference_assembly) + '.fai',
            ).base

        cram = alignment_input.resource_group(b)
        cram_file = cram[alignment_input.ext]

        # BAZAM requires indexed input.
        if not exists(alignment_input.index_path):
            index_cmd = f'samtools index {cram_file}'

        _reference_command_inp = (
            f'-Dsamjdk.reference_fasta={reference_inp}' if reference_inp else ''
        )

        bazam_cmd = dedent(
            f"""\
        bazam -Xmx16g {_reference_command_inp} \
        -n{min(nthreads, 6)} -bam {cram_file}{shard_param} \
        """
        )
        r1_param = 'r1'
        r2_param = ''
        fifo_commands[r1_param] = bazam_cmd

    else:
        assert isinstance(alignment_input, FastqPairs)
        use_bazam = False
        fastq_pairs = [p.as_resources(b) for p in alignment_input]
        files1 = [str(pair[0]) for pair in fastq_pairs]
        files2 = [str(pair[1]) for pair in fastq_pairs]
        if len(fastq_pairs) > 1:
            r1_param = 'r1'
            r2_param = 'r2'
            fifo_commands[r1_param] = f'cat {" ".join(files1)}'
            fifo_commands[r2_param] = f'cat {" ".join(files2)}'
        else:
            r1_param = files1[0]
            r2_param = files2[0]

    if aligner in [Aligner.BWAMEM2, Aligner.BWA]:
        if aligner == Aligner.BWAMEM2:
            tool_name = 'bwa-mem2'
            j.image(image_path('bwamem2'))
            index_exts = BWAMEM2_INDEX_EXTS
        else:
            tool_name = 'bwa'
            j.image(image_path('bwa'))
            index_exts = BWA_INDEX_EXTS
        bwa_reference = fasta_res_group(b, index_exts)
        rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
        # BWA command options:
        # -K   process INT input bases in each batch regardless of nThreads (for reproducibility)
        # -p   smart pairing (ignoring in2.fq)
        # -t16 threads
        # -Y   use soft clipping for supplementary alignments
        # -R   read group header line such as '@RG\tID:foo\tSM:bar'
        cmd = f"""\
        {tool_name} mem -K 100000000 \\
        {'-p' if use_bazam else ''} -t{nthreads - 1} -Y -R '{rg_line}' \\
        {bwa_reference.base} {r1_param} {r2_param}
        """

    elif aligner == Aligner.DRAGMAP:
        j.image(image_path('dragmap'))
        dragmap_index = b.read_input_group(
            **{
                k.replace('.', '_'): str(reference_path('broad/dragmap_prefix') / k)
                for k in DRAGMAP_INDEX_FILES
            }
        )
        if use_bazam:
            input_cmd = f'--interleaved=1 -b {r1_param}'
        else:
            input_cmd = f'-1 {r1_param} -2 {r2_param}'
        cmd = f"""\
        dragen-os -r {dragmap_index} {input_cmd} \\
            --RGID {sample_name} --RGSM {sample_name}
        """

    else:
        raise ValueError(f'Unsupported aligner: {aligner.value}')

    # prepare command for adding sort on the end
    cmd = dedent(cmd).strip()
    if should_sort:
        cmd += ' ' + sort_cmd(requested_nthreads) + f' -o {j.sorted_bam}'

    if fifo_commands:

        fifo_pre = [
            dedent(
                f"""
            mkfifo {fname}
            {command} > {fname} &
            pid_{fname}=$!
        """
            ).strip()
            for fname, command in fifo_commands.items()
        ]

        _fifo_waits = ' && '.join(
            f'wait $pid_{fname}' for fname in fifo_commands.keys()
        )
        fifo_post = dedent(
            f"""
            if {_fifo_waits}
            then
                echo -e "Background processes finished successfully"
            else
                # Background processes failed
                trap 'error1' ERR
            fi
        """
        ).strip()

        # Now prepare command
        cmd = '\n'.join([*fifo_pre, cmd, fifo_post])

    if index_cmd:
        cmd = dedent(index_cmd) + '\n' + cmd
    return j, cmd


def extract_fastq(
    b,
    cram: hb.ResourceGroup,
    ext: str = 'cram',
    job_attrs: dict | None = None,
    output_fq1: str | Path | None = None,
    output_fq2: str | Path | None = None,
) -> Job:
    """
    Job that converts a BAM or a CRAM file to an interleaved compressed fastq file.
    """
    j = b.new_job('Extract fastq', (job_attrs or {}) | dict(tool='bazam'))
    ncpu = 16
    nthreads = ncpu * 2  # multithreading
    j.cpu(ncpu)
    j.image(image_path('dragmap'))
    sequencing_type = get_config()['workflow']['sequencing_type']
    if sequencing_type == 'genome':
        j.storage('700G')

    reference = fasta_res_group(b)
    index_cmd = ''
    index_ext = 'crai' if ext == 'cram' else 'bai'
    if not index_ext not in cram:
        index_cmd = f'samtools index {cram[ext]}'
    cmd = f"""
    {index_cmd}
    bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
    -n{nthreads} -bam {cram[ext]} -r1 {j.fq1} -r2 {j.fq2}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if output_fq1 or output_fq2:
        assert output_fq1 and output_fq2, (output_fq1, output_fq2)
        b.write_output(j.fq1, str(output_fq1))
        b.write_output(j.fq2, str(output_fq2))
    return j


def sort_cmd(requested_nthreads: int) -> str:
    """
    Create command that coordinate-sorts SAM file
    """
    nthreads = STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()
    return dedent(
        f"""\
    | samtools sort -@{min(nthreads, 6) - 1} -T /io/batch/samtools-sort-tmp -Obam
    """
    ).strip()


def finalise_alignment(
    b: hb.Batch,
    align_cmd: str,
    stdout_is_sorted: bool,
    j: Job,
    requested_nthreads: int,
    markdup_tool: MarkDupTool,
    job_attrs: dict | None = None,
    output_path: CramPath | None = None,
    out_markdup_metrics_path: Path | None = None,
    align_cmd_out_fmt: str = 'sam',
    overwrite: bool = False,
) -> Job | None:
    """
    For `MarkDupTool.BIOBAMBAM`, adds bamsormadup command piped to the existing job.
    For `MarkDupTool.PICARD`, creates a new job, as Picard can't read from stdin.
    """

    reference = fasta_res_group(b)

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
        SO=coordinate M={j.markdup_metrics} outputformat=sam \\
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
            job_attrs=job_attrs,
            overwrite=overwrite,
        )

    if output_path:
        if md_j is not None:
            b.write_output(md_j.output_cram, str(output_path.path.with_suffix('')))
            if out_markdup_metrics_path:
                b.write_output(
                    md_j.markdup_metrics,
                    str(out_markdup_metrics_path),
                )
        else:
            b.write_output(j.sorted_bam, str(output_path.path.with_suffix('')))

    return md_j
