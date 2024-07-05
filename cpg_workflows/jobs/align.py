"""
FASTQ/BAM/CRAM -> CRAM: create Hail Batch jobs for (re-)alignment.
"""

import logging
import os.path
from ctypes import alignment
from enum import Enum
from math import e
from test.factories import dataset
from textwrap import dedent
from typing import cast
from weakref import ref

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import dataset_path, get_config, image_path, reference_path
from cpg_utils.hail_batch import command, fasta_res_group
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
    Get the reference used for the specific cram_version,
    so that Bazam is able to correctly decompress the reads.
    Note: Bazam is a tool that can output FASTQ in a form that
    can stream directly into common aligners such as BWA
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
    if get_config()['workflow'].get('realign_from_cram', False):
        logging.info('Realigning from CRAM')
        if (cram := sequencing_group.make_cram_path()).exists():
            logging.info(f'Realigning from CRAM {cram}')
            alignment_input = cram
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
    alignment_input: CramPath | BamPath,
    chr: str,
) -> Job:

    subset_cram_j = b.new_job('subset_tob_cram')
    subset_cram_j.image(image_path('samtools'))
    subset_cram_j.cpu(2)
    subset_cram_j.storage('150G')
    subset_cram_j.memory('standard')
    input_cram = b.read_input_group(**{'cram': alignment_input.path, 'crai': alignment_input.index_path})
    ref_path = fasta_res_group(b)['base']

    subset_cmd = f"""
    samtools view -T {ref_path} -C -o {subset_cram_j.output_cram} {input_cram['cram']} {chr} && \
    samtools index {subset_cram_j.output_cram} {subset_cram_j.output_crai}
    """
    subset_cram_j.output_cram.add_extension('.cram')
    subset_cram_j.output_crai.add_extension('.crai')
    subset_cram_j.command(command(subset_cmd))

    return subset_cram_j


def align(
    b,
    sequencing_group: SequencingGroup,
    job_attrs: dict | None = None,
    output_path: CramPath | None = None,
    out_markdup_metrics_path: Path | None = None,
    aligner: Aligner | str = Aligner.DRAGMAP,
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
        - is biobambam2, deduplication and alignment/merging are submitted within the same job.
        - is picard, deduplication is submitted in a separate job.

    - nthreads can be set for smaller test runs on toy instance, so the job
    doesn't take entire 32-cpu/64-threaded instance.
    """
    if isinstance(aligner, str):
        aligner = Aligner[aligner.upper()]

    if output_path and can_reuse(output_path.path, overwrite):
        return []

    # NOTE: If re-aligning from CRAM, the function call below returns a new CramPath
    # based on the [storage.<dataqset>] config key, ignoring the sequencing_group's
    # alignment input. The `index_path` attribute is set to `None` on this new instance
    # so `sharded_bazam` will be False. Not sure if this is a bug or intended behaviour?
    alignment_input = _get_alignment_input(sequencing_group)

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
    # sharded_bazam = (
    #     isinstance(alignment_input, CramPath | BamPath) and alignment_input.index_path and realignment_shards_num > 1
    # )
    sharded = sharded_fq  # or sharded_bazam

    jobs: list[Job] = []
    sharded_align_jobs = []
    sorted_bams = []

    if not sharded:  # Just running one alignment job
        if isinstance(alignment_input, FastqPairs):
            alignment_input = alignment_input[0]
        assert isinstance(alignment_input, FastqPair | BamPath | CramPath)
        logging.info(f'Aligning {alignment_input}. Either Fastq, Bam, or Cram. Not Sharded')
        if isinstance(alignment_input, BamPath | CramPath):
            assert alignment_input.index_path, alignment_input
            subset_cram_j = subset_cram(
                b,
                alignment_input,
                'chr21',
            )
            logging.info(f'Subsetted cram: {subset_cram_j}')
            alignment_input = CramPath(subset_cram_j.output_cram, subset_cram_j.output_crai).resource_group(b)
            assert isinstance(alignment_input, FastqPair | BamPath | CramPath)
        align_j, align_cmd = _align_one(
            b=b,
            job_name=base_job_name,
            alignment_input=alignment_input,
            requested_nthreads=requested_nthreads,
            sequencing_group_name=sequencing_group.id,
            extract_picard=get_config()['workflow'].get('extract_picard', False),
            job_attrs=job_attrs,
            aligner=aligner,
            should_sort=False,
        )
        stdout_is_sorted = False
        output_fmt = 'sam'
        jobs.append(align_j)
        logging.info(f'Alinging single job: {align_j}. No need to merge')
        merge_or_align_j = align_j

    else:  # Aligning in parallel and merging afterwards
        if sharded_fq:  # Aligning each lane separately, merging after
            logging.info(f'Aligning {alignment_input} fastq pairs. Sharded')
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

        else:
            assert isinstance(alignment_input, CramPath | BamPath)
            logging.info(f'Aligning {alignment_input} BAM or CRAM. Not Sharded')
            j, cmd = _align_one(
                b=b,
                job_name=base_job_name,
                alignment_input=alignment_input,
                sequencing_group_name=sequencing_group.id,
                job_attrs=job_attrs,
                aligner=aligner,
                requested_nthreads=requested_nthreads,
                should_sort=True,
            )
            # Sorting with samtools, but not adding deduplication yet, because we
            # need to merge first.
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
        sequencing_group=sequencing_group,
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
    b,
    job_name: str,
    alignment_input: FastqPair | CramPath | BamPath,
    requested_nthreads: int,
    sequencing_group_name: str,
    extract_picard: bool = False,
    job_attrs: dict | None = None,
    aligner: Aligner = Aligner.DRAGMAP,
    number_of_shards_for_realignment: int | None = None,
    shard_number: int | None = None,
    should_sort: bool = False,
) -> tuple[Job, str]:
    """
    Creates a job that (re)aligns reads to hg38. Returns the job object and a command
    separately, and doesn't add the command to the Job object, so stream-sorting
    and/or deduplication can be appended to the command later.

    Note: When this function is called within the align function, DRAGMAP is used as the default
    tool.
    """

    if number_of_shards_for_realignment is not None:
        assert number_of_shards_for_realignment > 1, number_of_shards_for_realignment

    job_attrs = (job_attrs or {}) | dict(label=job_name, tool=aligner.name)
    if shard_number is not None and number_of_shards_for_realignment is not None:
        job_name = f'{job_name} {shard_number + 1}/{number_of_shards_for_realignment} '
    job_name = f'{job_name} {alignment_input}'
    j = b.new_job(job_name, job_attrs)

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
    if isinstance(alignment_input, CramPath | BamPath):
        if extract_picard:
            logging.info("Using Picard to extract FASTQs from CRAM/BAM")
            extract_j = picard_extract_fastq(b, alignment_input.resource_group(b))
            fastq_pair = FastqPair(extract_j.fq1, extract_j.fq2).as_resources(b)
        else:
            # Extract fastqs from CRAM/BAM
            logging.info("Using samtools to extract FASTQs from CRAM/BAM")
            bam_or_cram_group = alignment_input.resource_group(b)
            extract_fastq_j = extract_fastq(b, bam_or_cram_group)
            fastq_pair = FastqPair(extract_fastq_j.fq1, extract_fastq_j.fq2).as_resources(b)
    else:  # only for BAMs that are missing index
        assert isinstance(alignment_input, FastqPair)
        fastq_pair = alignment_input.as_resources(b)

    r1_param = '$BATCH_TMPDIR/R1.fq.gz'
    r2_param = '$BATCH_TMPDIR/R2.fq.gz'
    # Need file names to end with ".gz" for BWA or DRAGMAP to parse correctly:
    prepare_fastq_cmd = dedent(
        f"""\
    mv {fastq_pair.r1} {r1_param}
    mv {fastq_pair.r2} {r2_param}
    """,
    )

    j.image(image_path('dragmap'))
    dragmap_index = b.read_input_group(
        **{k.replace('.', '_'): os.path.join(reference_path('broad/dragmap_prefix'), k) for k in DRAGMAP_INDEX_FILES},
    )
    input_params = f'-1 {r1_param} -2 {r2_param}'
    # TODO: consider reverting to use of all threads if node capacity
    # issue is resolved: https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/Job.20becomes.20unresponsive
    cmd = f"""\
    {prepare_fastq_cmd}
    dragen-os -r {dragmap_index} {input_params} \\
        --RGID {sequencing_group_name} --RGSM {sequencing_group_name} \\
        --num-threads {nthreads - 1}
    """

    if aligner is not Aligner.DRAGMAP:
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


def picard_extract_fastq(
    b: hb.Batch,
    bam_or_cram_group: hb.ResourceGroup,
    ext: str = 'cram',
    job_attrs: dict | None = None,
) -> Job:
    """
    Job that converts a BAM or a CRAM file to a pair of compressed fastq files.
    """
    collate_j = b.new_job('Collate BAM', (job_attrs or {}) | dict(tool='samtools'))
    collate_j.image(image_path('samtools'))
    extract_j = b.new_job('Extract fastq', (job_attrs or {}) | dict(tool='picard'))
    extract_j.image(image_path('picard'))
    reference_path = fasta_res_group(b)['base']
    res = STANDARD.request_resources(ncpu=16)
    if get_config()['workflow']['sequencing_type'] == 'genome':
        res.attach_disk_storage_gb = 700
    res.set_to_job(collate_j)
    res.set_to_job(extract_j)
    collate_j_cmd = f"""
    samtools collate --reference {reference_path} -@{res.get_nthreads() - 1} -u \
    -o {collate_j.collated_bam} {bam_or_cram_group[ext]}
    """
    collate_j.collated_bam.add_extension('.bam')
    collate_j.command(command(collate_j_cmd, monitor_space=True))
    extract_j_cmd = f"""
    picard SamToFastq I={collate_j.collated_bam} F=$BATCH_TMPDIR/R1.fq.gz F2=$BATCH_TMPDIR/R2.fq.gz
    mv $BATCH_TMPDIR/R1.fq.gz {extract_j.fq1}
    mv $BATCH_TMPDIR/R2.fq.gz {extract_j.fq2}
    """
    extract_j.command(command(extract_j_cmd, monitor_space=True))
    return extract_j


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
    # TODO: add ability to detect reference used in current BAM/CRAM and provide error handling
    reference_path = fasta_res_group(b)['base']
    res = STANDARD.request_resources(ncpu=16)
    if get_config()['workflow']['sequencing_type'] == 'genome':
        res.attach_disk_storage_gb = 700
    res.set_to_job(j)
    tmp_prefix = '$BATCH_TMPDIR/collate'
    cmd = f"""
    samtools collate --reference {reference_path} -@{res.get_nthreads() - 1} -u -O \
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
    """,
    ).strip()


def finalise_alignment(
    b: hb.Batch,
    sequencing_group: SequencingGroup,
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
    For `MarkDupTool.PICARD`, creates a new job, as Picard can't read from stdin.
    """

    reference = fasta_res_group(b)

    nthreads = STANDARD.request_resources(nthreads=requested_nthreads).get_nthreads()

    md_j = None
    align_cmd = align_cmd.strip()
    if not stdout_is_sorted:
        align_cmd += f' {sort_cmd(nthreads)}'
    align_cmd += f' > {j.sorted_bam}'

    j.command(command(align_cmd, monitor_space=True))  # type: ignore

    assert isinstance(j.sorted_bam, hb.ResourceFile)
    md_j = picard.markdup(
        b,
        j.sorted_bam,
        job_attrs=job_attrs,
        overwrite=overwrite,
    )

    if output_path:
        extract_picard = get_config()['workflow'].get('extract_picard', False)
        if extract_picard:
            output_path = CramPath(
                f'gs://cpg-tob-wgs-test/cram/picard_extracted_{sequencing_group.id}_chr21.cram',
            )
        else:
            output_path = CramPath(
                f'gs://cpg-tob-wgs-test/cram/samtools_extracted_{sequencing_group.id}_chr21.cram',
            )
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
