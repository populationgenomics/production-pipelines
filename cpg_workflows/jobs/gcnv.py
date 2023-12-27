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
    j = b.new_job(
        'Prepare intervals',
        job_attrs
        | {
            'tool': 'gatk PreprocessIntervals/AnnotateIntervals',
        },
    )
    j.image(image_path('gatk_gcnv'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(b)

    exclude_intervals = get_config()['workflow'].get('exclude_intervals', [])
    exclude_intervals_args = ' '.join(
        [f'--exclude-intervals {i}' for i in exclude_intervals]
    )

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
    j = b.new_job(
        'Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'}
    )
    j.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(b)

    j.declare_resource_group(
        counts={
            'counts.tsv.gz': '{root}.counts.tsv.gz',
            'counts.tsv.gz.tbi': '{root}.counts.tsv.gz.tbi',
        }
    )
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
        counts = b.read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            }
        )
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
    j = b.new_job(
        'Filter intervals and determine ploidy',
        job_attrs
        | {
            'tool': 'gatk FilterIntervals/DetermineGermlineContigPloidy',
        },
    )
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

        j = b.new_job(
            'Call germline CNVs',
            job_attrs
            | {
                'tool': 'gatk GermlineCNVCaller',
                'part': f'shard {i} of {n}',
            },
        )
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
    output_prefix: str,
) -> list[Job]:
    j = b.new_job(
        'Postprocess gCNV calls',
        job_attrs
        | {
            'tool': 'gatk PostprocessGermlineCNVCalls',
        },
    )
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
    allosomal_contigs_args = ' '.join(
        [f'--allosomal-contig {c}' for c in allosomal_contigs]
    )

    # declare all output files in advance
    j.declare_resource_group(
        output={
            'intervals.vcf.gz': '{root}/intervals.vcf.gz',
            'intervals.vcf.gz.tbi': '{root}/intervals.vcf.gz.tbi',
            'segments.vcf.gz': '{root}/segments.vcf.gz',
            'segments.vcf.gz.tbi': '{root}/segments.vcf.gz.tbi',
            'ratios.tsv': '{root}/ratios.tsv',
        }
    )

    postprocess_cmd = f"""
    gatk PostprocessGermlineCNVCalls \\
      --sequence-dictionary {reference.dict} {allosomal_contigs_args} \\
      --contig-ploidy-calls $BATCH_TMPDIR/ploidy-calls \\
      {model_shard_args} {calls_shard_args} \\
      --sample-index {sample_index} \\
      --output-genotyped-intervals {j.output['intervals.vcf.gz']} \\
      --output-genotyped-segments {j.output['segments.vcf.gz']} \\
      --output-denoised-copy-ratios {j.output['ratios.tsv']}
    """

    # index the output VCFs
    tabix_cmd = f"""
    tabix {j.output['intervals.vcf.gz']}
    tabix {j.output['segments.vcf.gz']}
    """

    j.command(command([*unpack_cmds, postprocess_cmd, tabix_cmd], setup_gcp=True))

    b.write_output(j.output, output_prefix)

    return [j]


def fix_intervals_vcf(
    b: hb.Batch, interval_vcf: Path, job_attrs: dict[str, str], output_path: Path
):
    """
    Note: the reheader loop is only required until the closure and
    adoption of https://github.com/broadinstitute/gatk/pull/8621
    Args:
        b (the batch instance):
        interval_vcf (Path): the individual intervals VCF
    Returns:
        the Job doing the work
    """
    reheader_job = b.new_job('Reheader intervals VCF', job_attrs | {'tool': 'bcftools'})
    reheader_job.declare_resource_group(
        output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    reheader_job.image(image_path('bcftools')).storage('1Gi')

    # read the Intervals VCF for this SG ID
    input_vcf = b.read_input(str(interval_vcf))

    # pull the header into a temp file
    reheader_job.command(f'bcftools view -h {input_vcf} > header')

    # sed command to swap Integer GT to String in-place
    reheader_job.command(
        r"sed -i 's/<ID=GT,Number=1,Type=Integer/<ID=GT,Number=1,Type=String/' header"
    )

    # apply the new header
    reheader_job.command(f'bcftools reheader -h header {input_vcf} -o temp.vcf.bgz')

    # split multiallelics (CNV calls are DEL/DUP at all loci)
    reheader_job.command(
        f'bcftools norm -m - temp.vcf.bgz | bgzip -c > {reheader_job.output["vcf.bgz"]}'
    )

    # and index with tabix
    reheader_job.command(f'tabix {reheader_job.output["vcf.bgz"]}')

    # get the output root to write to, and write both VCF and index
    b.write_output(reheader_job.output, str(output_path).removesuffix('.vcf.bgz'))

    return reheader_job


def merge_calls(
    b: hb.Batch,
    sg_vcfs: list[str],
    docker_image: str,
    job_attrs: dict[str, str],
    output_path: Path,
):
    """
    This job will run a fast simple merge on per-SGID call files
    It then throws in a python script to add in two additional header lines
    and edit the SVLEN and SVTYPE attributes into each row

    Args:
        b (batch):
        sg_vcfs (list[str]): paths to all individual VCFs
        docker_image (str): docker image to use
        job_attrs (dict): any params to atach to the job
        output_path (Path): path to the final merged VCF
    """

    if can_reuse(output_path):
        return None

    assert sg_vcfs, 'No VCFs to merge'

    j = b.new_job('Merge gCNV calls', job_attrs | {'tool': 'bcftools'})
    j.image(docker_image)

    # this should be made reactive, in case we scale past 10GB
    j.storage('10Gi')

    batch_vcfs = []
    for each_vcf in sg_vcfs:
        batch_vcfs.append(
            b.read_input_group(
                **{
                    'vcf.gz': each_vcf,
                    'vcf.gz.tbi': f'{each_vcf}.tbi',
                }
            )['vcf.gz']
        )

    j.declare_resource_group(
        output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: compression level
    j.command(
        f'bcftools merge {" ".join(batch_vcfs)} -Oz -o temp.vcf.bgz --threads 4 -m all -0'
    )
    j.command(
        fr"""
    python <<CODE
import gzip
headers = []
others = []
with gzip.open('temp.vcf.bgz', 'rt') as f:
    for line in f:
        if line.startswith('#'):
            headers.append(line)
            if line.startswith('##INFO=<ID=EN'):
                headers.extend([
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV Type">\n',
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV Length">\n']
                )
        else:
            l_split = line.split('\t')
            original_end = l_split[7]
            end_int = int(l_split[7].removeprefix('END='))
            alt_allele = l_split[4][1:-1]
            l_split[7] = 'SVTYPE={{alt}};SVLEN={{length}};{{end}}'.format(
                alt=alt_allele,
                length=str(end_int - int(l_split[1])),
                end=original_end
            )
            line = '\t'.join(l_split)
            others.append(line)
with open('temp.vcf', 'w') as f:
    f.writelines(headers)
    f.writelines(others)
CODE
    """
    )
    j.command(f'bgzip -c temp.vcf > {j.output["vcf.bgz"]}')
    j.command(f'tabix {j.output["vcf.bgz"]}')

    # get the output root to write to
    output_no_suffix = str(output_path).removesuffix('.vcf.bgz')
    b.write_output(j.output, output_no_suffix)
    return j
