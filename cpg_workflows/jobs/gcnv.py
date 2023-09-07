"""
Jobs that implement GATK-gCNV.
"""

from collections.abc import Iterable

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.resource import JobResourceFile, ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, fasta_res_group, image_path
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse


def prepare_intervals(
    b: hb.Batch,
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> list[Job]:
    j = b.new_job('Prepare intervals', job_attrs | {
        'tool': 'gatk PreprocessIntervals/AnnotateIntervals',
    })
    j.image(image_path('gatk_gcnv'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(b)

    exclude_intervals = get_config()['workflow'].get('exclude_intervals', [])
    exclude_intervals_args = ' '.join([f'--exclude-intervals {i}' for i in exclude_intervals])

    if sequencing_type == 'exome':
        intervals = b.read_input(get_config()['workflow'].get('intervals_path'))
        preprocess_cmd = f"""
        gatk PreprocessIntervals \\
          --reference {reference.base} --intervals {intervals} {exclude_intervals_args} \\
          --padding 250 --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY \\
          --output {j.preprocessed}
        """

    elif sequencing_type == 'genome':
        preprocess_cmd = f"""
        gatk PreprocessIntervals \\
          --reference {reference.base} {exclude_intervals_args} \\
          --interval-merging-rule OVERLAPPING_ONLY --padding 0 \\
          --output {j.preprocessed}
        """

    else:
        raise ValueError('workflow/sequencing_type must be one of: "exome", "genome"')

    annotate_cmd = f"""
    gatk AnnotateIntervals \\
      --reference {reference.base} --intervals {j.preprocessed} \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --output {j.annotated}
    """

    # Tell mypy the exact types of these Resources
    assert isinstance(j.preprocessed, JobResourceFile)
    assert isinstance(j.annotated, JobResourceFile)
    j.preprocessed.add_extension('.interval_list')
    j.annotated.add_extension('.tsv')

    j.command(command([preprocess_cmd, annotate_cmd]))
    for key, path in output_paths.items():
        b.write_output(j[key], str(path))
    return [j]


def collect_read_counts(
    b: hb.Batch,
    intervals_path: Path,
    cram_path: CramPath,
    job_attrs: dict[str, str],
    output_base_path: Path,
) -> list[Job]:
    j = b.new_job('Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'})
    j.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(b)

    j.declare_resource_group(counts={
        'counts.tsv.gz': '{root}.counts.tsv.gz',
        'counts.tsv.gz.tbi': '{root}.counts.tsv.gz.tbi',
    })
    assert isinstance(j.counts, ResourceGroup)

    cmd = f"""
    gatk CollectReadCounts \\
      --reference {reference.base} --intervals {intervals_path} \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --input {cram_path.path} --read-index {cram_path.index_path} \\
      --format TSV --output {j.counts}.counts.tsv

    bgzip {j.counts}.counts.tsv
    gatk IndexFeatureFile --input {j.counts['counts.tsv.gz']}
    """

    j.command(command(cmd, setup_gcp=True))
    b.write_output(j.counts, str(output_base_path))
    return [j]


def _counts_input_args(b: hb.Batch, counts_paths: Iterable[Path]) -> str:
    args = ''
    for f in counts_paths:
        counts = b.read_input_group(**{
            'counts.tsv.gz': str(f),
            'counts.tsv.gz.tbi': str(f) + '.tbi',
        })
        args += f' --input {counts["counts.tsv.gz"]}'

    return args


def filter_and_determine_ploidy(
    b: hb.Batch,
    ploidy_priors_path: Path,
    preprocessed_intervals_path: Path,
    annotated_intervals_path: Path,
    counts_paths: Iterable[Path],
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> list[Job]:
    j = b.new_job('Filter intervals and determine ploidy', job_attrs | {
        'tool': 'gatk FilterIntervals/DetermineGermlineContigPloidy',
    })
    j.image(image_path('gatk_gcnv'))

    counts_input_args = _counts_input_args(b, counts_paths)
    cmd = ''

    if can_reuse(output_paths['filtered']):
        # Remove 'filtered' entry from output_paths so we don't write_output() it later
        filtered: ResourceFile = b.read_input(str(output_paths.pop('filtered')))
    else:
        preprocessed_intervals = b.read_input(str(preprocessed_intervals_path))
        annotated_intervals = b.read_input(str(annotated_intervals_path))

        cmd += f"""
        gatk FilterIntervals \\
          --interval-merging-rule OVERLAPPING_ONLY \\
          --intervals {preprocessed_intervals} --annotated-intervals {annotated_intervals} \\
          {counts_input_args} \\
          --output {j.filtered}
        """

        assert isinstance(j.filtered, JobResourceFile)
        j.filtered.add_extension('.interval_list')
        filtered = j.filtered

    # (Other arguments may be cloud URLs, but this *must* be a local file)
    ploidy_priors = b.read_input(str(ploidy_priors_path))

    cmd += f"""
    gatk DetermineGermlineContigPloidy \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --intervals {filtered} --contig-ploidy-priors {ploidy_priors} \\
      {counts_input_args} \\
      --output $BATCH_TMPDIR --output-prefix ploidy

    tar -czf {j.calls} -C $BATCH_TMPDIR ploidy-calls
    tar -czf {j.model} -C $BATCH_TMPDIR ploidy-model
    """

    assert isinstance(j.calls, JobResourceFile)
    assert isinstance(j.model, JobResourceFile)
    j.calls.add_extension('.tar.gz')
    j.model.add_extension('.tar.gz')

    j.command(command(cmd))
    for key, path in output_paths.items():
        b.write_output(j[key], str(path))
    return [j]


def _shard_items(name_only: bool):
    if shard_partition := get_config()['workflow'].get('interval_shards'):
        count = len(shard_partition)
        for i, interval in enumerate(shard_partition, start=1):
            name = 'part' + 'c'.join([s.removeprefix('chr') for s in interval])
            if name_only:
                yield name
            else:
                chroms = '|'.join(interval)
                yield name, i, count, f"awk '/^@/ || $1 ~ /^({chroms})$/'"
    else:
        raise NotImplementedError('gCNV sharding technique not specified')


def shard_basenames():
    return _shard_items(name_only=True)


def shard_gcnv(
    b: hb.Batch,
    annotated_intervals_path: Path,
    filtered_intervals_path: Path,
    ploidy_calls_path: Path,
    counts_paths: Iterable[Path],
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> list[Job]:
    annotated_intervals = b.read_input(str(annotated_intervals_path))
    filtered_intervals = b.read_input(str(filtered_intervals_path))
    ploidy_calls_tarball = b.read_input(str(ploidy_calls_path))
    counts_input_args = _counts_input_args(b, counts_paths)

    jobs: list[Job] = []

    for name, i, n, select_cmd in _shard_items(name_only=False):
        if can_reuse(output_paths[name]):
            continue

        j = b.new_job('Call germline CNVs', job_attrs | {
            'tool': 'gatk GermlineCNVCaller',
            'part': f'shard {i} of {n}',
        })
        j.image(image_path('gatk_gcnv'))
        j.memory('16Gi')  # TODO revisit limits

        cmd = f"""
        tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR

        {select_cmd} < {filtered_intervals} > {j.shard_intervals}

        gatk GermlineCNVCaller \\
          --run-mode COHORT --interval-merging-rule OVERLAPPING_ONLY \\
          --intervals {j.shard_intervals} --annotated-intervals {annotated_intervals} \\
          {counts_input_args} \\
          --contig-ploidy-calls $BATCH_TMPDIR/ploidy-calls \\
          --output $BATCH_TMPDIR --output-prefix {name}

        tar -czf {j.shard_tarball} -C $BATCH_TMPDIR {name}-calls {name}-model {name}-tracking
        """

        assert isinstance(j.shard_intervals, JobResourceFile)
        assert isinstance(j.shard_tarball, JobResourceFile)
        j.shard_intervals.add_extension('.interval_list')
        j.shard_tarball.add_extension('.tar.gz')

        j.command(command(cmd))
        b.write_output(j.shard_tarball, str(output_paths[name]))
        jobs.append(j)

    return jobs


def postprocess_calls(
    b: hb.Batch,
    ploidy_calls_path: Path,
    shard_paths: dict[str, Path],
    sample_index: int,
    job_attrs: dict[str, str],
    output_path: dict[str, Path],
) -> list[Job]:
    j = b.new_job('Postprocess gCNV calls', job_attrs | {
        'tool': 'gatk PostprocessGermlineCNVCalls',
    })
    j.image(image_path('gatk_gcnv'))
    j.storage('12Gi')  # TODO revisit limits

    reference = fasta_res_group(b)

    ploidy_calls_tarball = b.read_input(str(ploidy_calls_path))

    unpack_cmds = [f'tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR']

    model_shard_args = ''
    calls_shard_args = ''
    for name, path in shard_paths.items():
        unpack_cmds.append(f'gsutil cat {path} | tar -xz -C $BATCH_TMPDIR')
        model_shard_args += f' --model-shard-path $BATCH_TMPDIR/{name}-model'
        calls_shard_args += f' --calls-shard-path $BATCH_TMPDIR/{name}-calls'

    allosomal_contigs = get_config()['workflow'].get('allosomal_contigs', [])
    allosomal_contigs_args = ' '.join([f'--allosomal-contig {c}' for c in allosomal_contigs])

    postprocess_cmd = f"""
    gatk PostprocessGermlineCNVCalls \\
      --sequence-dictionary {reference.dict} {allosomal_contigs_args} \\
      --contig-ploidy-calls $BATCH_TMPDIR/ploidy-calls \\
      {model_shard_args} {calls_shard_args} \\
      --sample-index {sample_index} \\
      --output-genotyped-intervals {j.intervals} \\
      --output-genotyped-segments {j.segments} \\
      --output-denoised-copy-ratios {j.ratios}
    """

    assert isinstance(j.intervals, JobResourceFile)
    assert isinstance(j.segments, JobResourceFile)
    assert isinstance(j.ratios, JobResourceFile)
    j.intervals.add_extension('.vcf.gz')
    j.segments.add_extension('.vcf.gz')
    j.ratios.add_extension('.tsv')

    j.command(command([*unpack_cmds, postprocess_cmd], setup_gcp=True))
    for key, path in output_path.items():
        b.write_output(j[key], str(path))
    return [j]
