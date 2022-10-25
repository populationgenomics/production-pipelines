"""
FASTQ/BAM/CRAM -> CRAM: create Hail Batch jobs for (re-)alignment.
"""

from enum import Enum
from textwrap import dedent
from typing import cast
import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from cpg_utils.hail_batch import command
from cpg_utils.workflows.targets import Sample
from cpg_utils.workflows.filetypes import (
    AlignmentInput,
    FastqPairs,
    CramPath,
    BamPath,
    FastqPair,
)
from cpg_utils.workflows.resources import STANDARD
from cpg_utils.workflows.utils import can_reuse

from . import picard


BWA_INDEX_EXTS = ['sa', 'amb', 'bwt', 'ann', 'pac', 'alt']
BWAMEM2_INDEX_EXTS = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac', 'alt']
DRAGMAP_INDEX_FILES = ['hash_table.cfg.bin', 'hash_table.cmp', 'reference.bin']


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
    so that Bazam is able to correctly decompress the reads
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
            logging.info(f'Realigning from {realign_cram_ver} CRAM {path}')
            alignment_input = CramPath(
                path,
                reference_assembly=_get_cram_reference_from_version(realign_cram_ver),
            )

    if not alignment_input:
        raise MissingAlignmentInputException(
            f'No alignment inputs found for sample {sample}'
            + (f': {alignment_input}' if alignment_input else '')
        )

    if get_config()['workflow'].get('check_inputs', True):
        if not alignment_input.exists():
            raise MissingAlignmentInputException(
                f'Alignment inputs for sample {sample} do not exist '
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

    base_job_name = 'Align'
    if extra_label:
        base_job_name += f' {extra_label}'

    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_threads()

    if get_config()['workflow']['sequencing_type'] == 'genome':
        realignment_shards_num = 10
    else:
        assert get_config()['workflow']['sequencing_type'] == 'exome'
        realignment_shards_num = 1

    sharded_fq = isinstance(alignment_input, FastqPairs) and len(alignment_input) > 1
    sharded_bazam = (
        isinstance(alignment_input, CramPath | BamPath)
        and alignment_input.index_path
        and realignment_shards_num > 1
    )
    sharded = sharded_fq or sharded_bazam

    jobs: list[Job] = []
    sharded_align_jobs = []
    sorted_bams = []

    if not sharded:  # Just running one alignment job
        if isinstance(alignment_input, FastqPairs):
            alignment_input = alignment_input[0]
        align_j, align_cmd = _align_one(
            b=b,
            job_name=base_job_name,
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
                    job_name=base_job_name,
                    alignment_input=pair,
                    requested_nthreads=requested_nthreads,
                    sample_name=sample.id,
                    job_attrs=job_attrs,
                    aligner=aligner,
                    should_sort=True,
                )
                assert isinstance(j, Job)
                j.command(command(cmd, monitor_space=True))  # type: ignore
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        elif sharded_bazam:  # Using BAZAM to shard CRAM
            assert realignment_shards_num, realignment_shards_num
            for shard_number in range(realignment_shards_num):
                j, cmd = _align_one(
                    b=b,
                    job_name=base_job_name,
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
                j.command(command(cmd, monitor_space=True))  # type: ignore
                sorted_bams.append(j.sorted_bam)
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
        samtools merge -@{nthreads - 1} - {' '.join(map(str, sorted_bams))}
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
    if md_j and md_j != merge_or_align_j:
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
        if isinstance(alignment_input, FastqPair | FastqPairs | BamPath):
            storage_gb = 400  # for WGS FASTQ or BAM inputs, more disk is needed
    return storage_gb


def _align_one(
    b,
    job_name: str,
    alignment_input: FastqPair | CramPath | BamPath,
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
    job_name = f'{job_name} {alignment_input}'
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

    if isinstance(alignment_input, CramPath | BamPath):
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
        if isinstance(alignment_input, CramPath):
            assert (
                alignment_input.reference_assembly
            ), f'The reference input for the alignment input "{alignment_input.path}" was not set'
            reference_inp = b.read_input_group(
                base=str(alignment_input.reference_assembly),
                fai=str(alignment_input.reference_assembly) + '.fai',
            ).base

        group = alignment_input.resource_group(b)

        _reference_command_inp = (
            f'-Dsamjdk.reference_fasta={reference_inp}' if reference_inp else ''
        )

        sort_index_input_cmd = ''
        if not alignment_input.index_path:
            sort_index_input_cmd = dedent(
                f"""
            mkdir -p $BATCH_TMPDIR/sorted
            mkdir -p $BATCH_TMPDIR/sort_tmp
            samtools sort \
            {group[alignment_input.ext]} \
            -@{nthreads - 1} \
            -T $BATCH_TMPDIR/sort_tmp \
            > $BATCH_TMPDIR/sorted.{alignment_input.ext}
            
            mv $BATCH_TMPDIR/sorted.{alignment_input.ext} {group[alignment_input.ext]}
            rm -rf $BATCH_TMPDIR/sort_tmp
            
            samtools index -@{nthreads - 1} {group[alignment_input.ext]}
            """
            )

        bazam_cmd = dedent(
            f"""\
        {sort_index_input_cmd}
        bazam -Xmx16g {_reference_command_inp} \
        -n{min(nthreads, 6)} -bam {group[alignment_input.ext]}{shard_param} \
        """
        )
        prepare_fastq_cmd = ''
        r1_param = 'r1'
        r2_param = ''
        fifo_commands[r1_param] = bazam_cmd

    else:  # only for BAMs that are missing index
        assert isinstance(alignment_input, FastqPair)
        fastq_pair = alignment_input.as_resources(b)
        use_bazam = False
        r1_param = '$BATCH_TMPDIR/R1.fq.gz'
        r2_param = '$BATCH_TMPDIR/R2.fq.gz'
        # Need file names to end with ".gz" for BWA or DRAGMAP to parse correctly:
        prepare_fastq_cmd = dedent(
            f"""\
        mv {fastq_pair.r1} {r1_param}
        mv {fastq_pair.r2} {r2_param}
        """
        )

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
        {prepare_fastq_cmd}
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
            input_params = f'--interleaved=1 -b {r1_param}'
        else:
            input_params = f'-1 {r1_param} -2 {r2_param}'
        cmd = f"""\
        {prepare_fastq_cmd}
        dragen-os -r {dragmap_index} {input_params} \\
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
            {cmd} > {fname} &
            pid_{fname}=$!
        """
            ).strip()
            for fname, cmd in fifo_commands.items()
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
    return j, cmd


def extract_fastq(
    b,
    bam_or_cram_group: hb.ResourceGroup,
    ext: str = 'cram',
    job_attrs: dict | None = None,
    output_fq1: str | Path | None = None,
    output_fq2: str | Path | None = None,
) -> Job:
    """
    Job that converts a BAM or a CRAM file to a pair of compressed fastq files.
    """
    j = b.new_job('Extract fastq', (job_attrs or {}) | dict(tool='samtools'))
    j.image(image_path('samtools'))
    res = STANDARD.request_resources(ncpu=16)
    if get_config()['workflow']['sequencing_type'] == 'genome':
        res.attach_disk_storage_gb = 700
    res.set_to_job(j)
    tmp_prefix = '$BATCH_TMPDIR/collate'
    cmd = f"""
    samtools collate -@{res.get_nthreads() - 1} -u -O \
    {bam_or_cram_group[ext]} {tmp_prefix} | \\
    samtools fastq -@{res.get_nthreads() - 1} \
    -1 $BATCH_TMPDIR/R1.fq.gz -2 $BATCH_TMPDIR/R2.fq.gz \
    -0 /dev/null -s /dev/null -n
    # Can't write directly to j.fq1 and j.fq2 because samtools-fastq requires the 
    # file names to end with ".gz" in order to create compressed outputs.
    mv $BATCH_TMPDIR/R1.fq.gz {j.fq1}
    mv $BATCH_TMPDIR/R2.fq.gz {j.fq2}
    """
    j.command(command(cmd, monitor_space=True))
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
    | samtools sort -@{min(nthreads, 6) - 1} -T $BATCH_TMPDIR/samtools-sort-tmp -Obam
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
        j.declare_resource_group(  # type: ignore
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            }
        )
        assert isinstance(j.output_cram, hb.ResourceGroup)
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

    j.command(command(align_cmd, monitor_space=True))  # type: ignore

    assert isinstance(j.sorted_bam, hb.ResourceFile)
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
