"""
FASTQ/BAM/CRAM -> CRAM: create Hail Batch jobs for (re-)alignment.
"""

import logging
import os.path
from enum import Enum
from textwrap import dedent
from typing import cast

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import command, fasta_res_group, output_path
from cpg_workflows.filetypes import (
    AlignmentInput,
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_workflows.resources import HIGHMEM, STANDARD
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import can_reuse

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
    Get the reference used for the specific cram_version
    """
    cram_version_map = get_config()['workflow'].get('cram_version_reference', {})
    if cram_version in cram_version_map:
        return cram_version_map[cram_version]
    raise ValueError(
        f'Unrecognised cram_version: "{cram_version}", expected one of: {", ".join(cram_version_map.keys())}',
    )


class MissingAlignmentInputException(Exception):
    """Raise if alignment input is missing"""

    pass


def _get_alignment_input(sequencing_group: SequencingGroup) -> AlignmentInput:
    """Given a sequencing group, will return an AlignmentInput object that
    represents the path to a relevant input (e.g. CRAM/BAM path)"""
    sequencing_type = get_config()['workflow']['sequencing_type']
    alignment_input = sequencing_group.alignment_input_by_seq_type.get(sequencing_type)
    if realign_cram_ver := get_config()['workflow'].get('realign_from_cram_version'):
        if (
            path := (sequencing_group.dataset.prefix() / 'cram' / realign_cram_ver / f'{sequencing_group.id}.cram')
        ).exists():
            logging.info(f'Realigning from {realign_cram_ver} CRAM {path}')
            alignment_input = CramPath(
                path,
                reference_assembly=_get_cram_reference_from_version(realign_cram_ver),
            )
    if realign_from_cram := get_config()['workflow'].get('realign_from_cram'):
        alignment_input = CramPath(
            path=(sequencing_group.dataset.prefix() / 'cram' / f'{sequencing_group.id}.cram'),
            reference_assembly=get_config()['workflow']['ref_fasta'],
        )
        logging.info(f'Realigning from CRAM {alignment_input.path}')
    if not alignment_input:
        raise MissingAlignmentInputException(
            f'No alignment inputs found for sequencing group {sequencing_group}'
            + (f': {alignment_input}' if alignment_input else ''),
        )

    if get_config()['workflow'].get('check_inputs', True):
        if not alignment_input.exists():
            raise MissingAlignmentInputException(
                f'Alignment inputs for sequencing group {sequencing_group} do not exist '
                + (f': {alignment_input}' if alignment_input else ''),
            )

    return alignment_input


def subset_cram(
    b: hb.Batch,
    bam_or_cram_group: hb.ResourceGroup,
    alignment_input: FastqPair | CramPath | BamPath,
    sequencing_group: SequencingGroup,
    chr: str,
    output_bucket: str = 'tmp',
) -> Job:
    logging.info(f'alignment_input: {type(alignment_input)}')
    subset_cram_j = b.new_job('subset_tob_cram')
    subset_cram_j.image(image_path('samtools'))
    subset_cram_j.cpu(2)
    subset_cram_j.storage('150G')
    subset_cram_j.memory('standard')
    if isinstance(alignment_input, FastqPair):
        ref_path = fasta_res_group(b)['base']
        logging.info(f'Reference path: {ref_path}')
    if isinstance(alignment_input, CramPath | BamPath):
        ref_path = b.read_input(alignment_input.reference_assembly)
        logging.info(f'Reference path: {ref_path}')
    subset_cram_j.declare_resource_group(
        cram_output={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        },
    )
    subset_cmd = f"""
    samtools view -T {ref_path} -C -o {subset_cram_j.cram_output.cram} {bam_or_cram_group['cram']} {chr} && \
    samtools index {subset_cram_j.cram_output.cram}
    """
    subset_cram_j.command(command(subset_cmd))

    if output_bucket:
        path = f'subset_{chr}/{sequencing_group.id}_{chr}'
        logging.info(f'Writing subset cram output to {output_path(path, "tmp")}')
        b.write_output(subset_cram_j.cram_output, output_path(path, 'tmp'))

    return subset_cram_j


def align(
    b,
    sequencing_group: SequencingGroup,
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

        - for dragmap, submit an extra job to extract a pair of fastqs from the cram/bam

    - if the markdup tool:
        - is biobambam2, deduplication and alignment/merging are submitted within the same job.
        - is picard, deduplication is submitted in a separate job.

    - nthreads can be set for smaller test runs on toy instance, so the job
    doesn't take entire 32-cpu/64-threaded instance.
    """
    if output_path and can_reuse(output_path.path, overwrite):
        return []

    # NOTE: If re-aligning from CRAM, the function call below returns a new CramPath
    # based on the [storage.<dataqset>] config key, ignoring the sequencing_group's
    # alignment input. The `index_path` attribute is set to `None` on this new instance
    alignment_input = _get_alignment_input(sequencing_group)

    base_job_name = 'Align'
    if extra_label:
        base_job_name += f' {extra_label}'

    # if number of threads is not requested, using whole instance
    requested_nthreads = requested_nthreads or STANDARD.max_threads()

    sharded_fq = isinstance(alignment_input, FastqPairs) and len(alignment_input) > 1

    jobs: list[Job] = []
    sharded_align_jobs = []
    sorted_bams = []

    if not sharded_fq:  # Just running one alignment job
        if isinstance(alignment_input, FastqPairs):
            alignment_input = alignment_input[0]
        assert isinstance(alignment_input, FastqPair | BamPath | CramPath)
        extract_reads = False
        if isinstance(alignment_input, CramPath | BamPath):
            extract_reads = True
            if isinstance(alignment_input, CramPath):
                assert (
                    alignment_input.reference_assembly
                ), f'The reference input for the alignment input "{alignment_input.path}" was not set'
        align_j, align_cmd = _align_one(
            b=b,
            job_name=base_job_name,
            alignment_input=alignment_input,  # type: ignore[arg-type]
            requested_nthreads=requested_nthreads,
            sequencing_group_name=sequencing_group.id,
            job_attrs=job_attrs,
            aligner=aligner,
            should_sort=False,
            extract_reads=extract_reads,
        )
        stdout_is_sorted = False
        output_fmt = 'sam'
        jobs.append(align_j)
        logging.info('Aligning single job. No need to merge')
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
                    sequencing_group_name=sequencing_group.id,
                    job_attrs=job_attrs,
                    aligner=aligner,
                    should_sort=True,
                )
                assert isinstance(j, Job)
                j.command(command(cmd, monitor_space=True))  # type: ignore
                sorted_bams.append(j.sorted_bam)
                sharded_align_jobs.append(j)

        merge_j = b.new_job('Merge BAMs', (job_attrs or {}) | dict(tool='samtools_merge'))
        merge_j.image(image_path('samtools'))

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
    This function pulls the value from the config, otherwise utilises the full capacity of a standard
    machine. Increases storage if genomes are being handled, and again for unindexed/unsorted CRAM/BAMs.
    """
    storage_gb = None  # avoid attaching extra disk by default

    try:
        storage_gb = get_config()['workflow']['resources']['Align']['storage_gb']
    except KeyError:
        pass
    else:
        assert isinstance(storage_gb, int), storage_gb
        return storage_gb

    # Taking a full instance without attached by default:
    storage_gb = STANDARD.calc_instance_disk_gb()
    if get_config()['workflow']['sequencing_type'] == 'genome':
        # More disk is needed for FASTQ or BAM inputs than for realignment from CRAM
        if isinstance(alignment_input, FastqPair | FastqPairs | BamPath):
            storage_gb = 400
        # For unindexed/unsorted CRAM or BAM inputs, extra storage is needed for tmp
        if isinstance(alignment_input, CramPath | BamPath) and not alignment_input.index_path:
            storage_gb += 150
    return storage_gb


def _align_one(
    b: hb.Batch,
    job_name: str,
    alignment_input: FastqPair | CramPath | BamPath,
    requested_nthreads: int,
    sequencing_group_name: str,
    job_attrs: dict | None = None,
    aligner: Aligner = Aligner.BWA,
    should_sort: bool = False,
    extract_reads: bool = False,
) -> tuple[Job, str]:
    """
    Creates a job that (re)aligns reads to hg38. Returns the job object and a command
    separately, and doesn't add the command to the Job object, so stream-sorting
    and/or deduplication can be appended to the command later.

    Note: When this function is called within the align function, DRAGMAP is used as the default
    tool.
    """
    job_attrs = (job_attrs or {}) | dict(label=job_name, tool=aligner.name)
    job_name = f'{job_name} {alignment_input}'
    j: Job = b.new_job(job_name, job_attrs)

    if get_config()['resource_overrides'].get('align_use_highmem'):
        align_machine_type = HIGHMEM
    else:
        align_machine_type = STANDARD

    nthreads = align_machine_type.set_resources(
        j,
        nthreads=requested_nthreads,
        storage_gb=storage_for_align_job(alignment_input=alignment_input),
    ).get_nthreads()

    sort_index_input_cmd = ''

    # 2022-07-22 mfranklin:
    #   Replace process substitution with named-pipes (FIFO)
    #   This is named-pipe name -> command to populate it
    fifo_commands: dict[str, str] = {}
    use_interleaved = False
    if isinstance(alignment_input, CramPath | BamPath):
        if extract_reads:
            logging.info("Using samtools to extract FASTQs from CRAM/BAM")
            logging.info(f"Alignment input: {alignment_input.path} {alignment_input.index_path}")
            bam_or_cram_group = alignment_input.resource_group(b)
            extract_fastq_j = extract_fastq(
                b,
                bam_or_cram_group,
                alignment_input,
                ext=alignment_input.ext,
            )
            j.depends_on(extract_fastq_j)
            interleaved_fastq = extract_fastq_j.all_reads
            interleave_param = '$BATCH_TMPDIR/all_reads.fq.gz'
            # Need file names to end with ".gz" for DRAGMAP to parse correctly:
            prepare_fastq_cmd = dedent(
                f"""\
            mv {interleaved_fastq} {interleave_param}
            """,
            )
            use_interleaved = True

        prepare_fastq_cmd = ''
        r1_param = 'r1'
        r2_param = ''

    else:  # only for BAMs that are missing index
        assert isinstance(alignment_input, FastqPair)
        fastq_pair = alignment_input.as_resources(b)
        r1_param = '$BATCH_TMPDIR/R1.fq.gz'
        r2_param = '$BATCH_TMPDIR/R2.fq.gz'
        # Need file names to end with ".gz" for DRAGMAP to parse correctly:
        prepare_fastq_cmd = dedent(
            f"""\
        mv {fastq_pair.r1} {r1_param}
        mv {fastq_pair.r2} {r2_param}
        """,
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
        rg_line = f'@RG\\tID:{sequencing_group_name}\\tSM:{sequencing_group_name}'
        # BWA command options:
        # -K   process INT input bases in each batch regardless of nThreads (for reproducibility)
        # -p   smart pairing (ignoring in2.fq)
        # -t16 threads
        # -Y   use soft clipping for supplementary alignments
        # -R   read group header line such as '@RG\tID:foo\tSM:bar'
        cmd = f"""\
        {prepare_fastq_cmd}
        {tool_name} mem -K 100000000 \\
        -t{nthreads - 1} -Y -R '{rg_line}' \\
        {bwa_reference.base} {r1_param} {r2_param}
        """

    elif aligner == Aligner.DRAGMAP:
        j.image(image_path('dragmap'))
        dragmap_index = b.read_input_group(
            **{
                k.replace('.', '_'): os.path.join(reference_path('broad/dragmap_prefix'), k)
                for k in DRAGMAP_INDEX_FILES
            },
        )
        if use_interleaved:
            input_params = f'--interleaved=1 -b {interleave_param}'
        else:
            input_params = f'-1 {r1_param} -2 {r2_param}'
        # TODO: consider reverting to use of all threads if node capacity
        # issue is resolved: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/Job.20becomes.20unresponsive
        cmd = f"""\
        {prepare_fastq_cmd}
        dragen-os -r {dragmap_index} {input_params} \\
            --RGID {sequencing_group_name} --RGSM {sequencing_group_name} \\
            --num-threads {nthreads - 1}
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
        """,
            ).strip()
            for fname, cmd in fifo_commands.items()
        ]

        _fifo_waits = ' && '.join(f'wait $pid_{fname}' for fname in fifo_commands.keys())
        fifo_post = dedent(
            f"""
            if {_fifo_waits}
            then
                echo -e "Background processes finished successfully"
            else
                # Background processes failed
                trap 'error1' ERR
            fi
        """,
        ).strip()

        # Now prepare command
        cmd = '\n'.join([sort_index_input_cmd, *fifo_pre, cmd, fifo_post])
    return j, cmd


def extract_fastq(
    b,
    bam_or_cram_group: hb.ResourceGroup,
    alignment_input: CramPath | BamPath,
    ext: str = 'cram',
    job_attrs: dict | None = None,
) -> Job:
    """
    Job that converts a BAM or a CRAM file to a pair of compressed fastq files.
    """
    j = b.new_job('Extract fastq', (job_attrs or {}) | dict(tool='samtools'))
    j.image(image_path('samtools'))
    # TODO: add ability to detect reference used in current BAM/CRAM and provide error handling
    if isinstance(alignment_input, CramPath):
        reference_flag = f'--reference f{b.read_input(str(alignment_input.reference_assembly))}'
    else:
        reference_flag = ''
    res = STANDARD.request_resources(ncpu=16)
    if get_config()['workflow']['sequencing_type'] == 'genome':
        res.attach_disk_storage_gb = 700
    res.set_to_job(j)
    tmp_prefix = '$BATCH_TMPDIR/collate'
    cmd = f"""
    samtools collate {reference_flag} -@{res.get_nthreads() - 1} -u -O \
    {bam_or_cram_group[ext]} {tmp_prefix} | \\
    samtools fastq -@{res.get_nthreads() - 1} \
    -0 /dev/null -n | gzip > $BATCH_TMPDIR/all_reads.fq.gz
    # Can't write directly to j.all_reads because samtools-fastq requires the
    # file names to end with ".gz" in order to create compressed outputs.
    mv $BATCH_TMPDIR/all_reads.fq.gz {j.all_reads}
    """
    j.command(command(cmd, monitor_space=True))

    return j


def sort_cmd(requested_nthreads: int) -> str:
    """
    Create command that coordinate-sorts SAM file
    """
    nthreads = STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()
    return dedent(
        f"""\
    | samtools sort -@{min(nthreads, 6) - 1} -T $BATCH_TMPDIR/samtools-sort-tmp -Obam
    """,
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
            },
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
