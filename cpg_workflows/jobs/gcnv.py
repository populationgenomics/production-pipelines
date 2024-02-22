"""
Jobs that implement GATK-gCNV.
"""

from collections.abc import Iterable

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.resource import (
    JobResourceFile,
    ResourceFile,
    ResourceGroup,
    Resource,
)

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, fasta_res_group, get_batch, image_path
from cpg_workflows.filetypes import CramPath
from cpg_workflows.stages.gatk_sv.gatk_sv_common import reference_path
from cpg_workflows.utils import can_reuse, chunks


def prepare_intervals(
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> Job:
    j = get_batch().new_job(
        'Prepare intervals',
        job_attrs
        | {
            'tool': 'gatk PreprocessIntervals/AnnotateIntervals',
        },
    )
    j.image(image_path('gatk_gcnv'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(get_batch())

    exclude_intervals = get_config()['workflow'].get('exclude_intervals', [])
    exclude_intervals_args = ' '.join(
        [f'--exclude-intervals {i}' for i in exclude_intervals]
    )

    if sequencing_type == 'exome':
        intervals = get_batch().read_input(
            str(reference_path('broad/exome_calling_interval_lists'))
        )
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
        get_batch().write_output(j[key], str(path))
    return j


def collect_read_counts(
    intervals_path: Path,
    cram_path: CramPath,
    job_attrs: dict[str, str],
    output_base_path: Path,
) -> list[Job]:
    j = get_batch().new_job(
        'Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'}
    )
    j.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(get_batch())

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
    get_batch().write_output(j.counts, str(output_base_path))
    return [j]


def _counts_input_args(counts_paths: Iterable[Path]) -> str:
    args = ''
    for f in counts_paths:
        counts = get_batch().read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            }
        )
        args += f' --input {counts["counts.tsv.gz"]}'

    return args


def filter_and_determine_ploidy(
    ploidy_priors_path: str,
    preprocessed_intervals_path: Path,
    annotated_intervals_path: Path,
    counts_paths: Iterable[Path],
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> list[Job]:
    j = get_batch().new_job(
        'Filter intervals and determine ploidy',
        job_attrs
        | {
            'tool': 'gatk FilterIntervals/DetermineGermlineContigPloidy',
        },
    )
    j.image(image_path('gatk_gcnv'))

    counts_input_args = _counts_input_args(counts_paths)
    cmd = ''

    if can_reuse(output_paths['filtered']):
        # Remove 'filtered' entry from output_paths so we don't write_output() it later
        filtered: ResourceFile = get_batch().read_input(
            str(output_paths.pop('filtered'))
        )
    else:
        preprocessed_intervals = get_batch().read_input(
            str(preprocessed_intervals_path)
        )
        annotated_intervals = get_batch().read_input(str(annotated_intervals_path))

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
    ploidy_priors = get_batch().read_input(ploidy_priors_path)

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
        get_batch().write_output(j[key], str(path))
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
    annotated_intervals_path: Path,
    filtered_intervals_path: Path,
    ploidy_calls_path: Path,
    counts_paths: Iterable[Path],
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> list[Job]:
    annotated_intervals = get_batch().read_input(str(annotated_intervals_path))
    filtered_intervals = get_batch().read_input(str(filtered_intervals_path))
    ploidy_calls_tarball = get_batch().read_input(str(ploidy_calls_path))
    counts_input_args = _counts_input_args(counts_paths)

    jobs: list[Job] = []

    for name, i, n, select_cmd in _shard_items(name_only=False):
        if can_reuse(output_paths[name]):
            continue

        j = get_batch().new_job(
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
        get_batch().write_output(j.shard_tarball, str(output_paths[name]))
        jobs.append(j)

    return jobs


def postprocess_calls_1(
    ploidy_calls_path: Path,
    shard_paths: dict[str, Path],
    sample_index: int,
    job_attrs: dict[str, str],
    output_prefix: str,
) -> Job:
    j = get_batch().new_job(
        'Postprocess gCNV calls',
        job_attrs
        | {
            'tool': 'gatk PostprocessGermlineCNVCalls',
        },
    )
    j.image(image_path('gatk_gcnv'))
    j.storage('12Gi')  # TODO revisit limits

    reference = fasta_res_group(get_batch())

    ploidy_calls_tarball = get_batch().read_input(str(ploidy_calls_path))

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

    get_batch().write_output(j.output, output_prefix)

    return j


def fix_intervals_vcf(interval_vcf: Path, job_attrs: dict[str, str], output_path: Path):
    """
    Note: the reheader loop is only required until the closure and
    adoption of https://github.com/broadinstitute/gatk/pull/8621
    Args:
        interval_vcf (Path): the individual intervals VCF
    Returns:
        the Job doing the work
    """
    reheader_job = get_batch().new_job(
        'Reheader intervals VCF', job_attrs | {'tool': 'bcftools'}
    )
    reheader_job.declare_resource_group(
        output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    reheader_job.image(image_path('bcftools')).storage('1Gi')

    # read the Intervals VCF for this SG ID
    input_vcf = get_batch().read_input(str(interval_vcf))

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
    get_batch().write_output(
        reheader_job.output, str(output_path).removesuffix('.vcf.bgz')
    )

    return reheader_job


def joint_segment_vcfs(
    segment_vcfs: list[ResourceFile],
    pedigree: ResourceFile,
    reference: ResourceGroup,
    intervals: ResourceFile,  # hmm,
    job_attrs: dict,
) -> tuple[Job, Resource]:
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates

    Returns:
        the job that does the work, and the resulting resource group of VCF & index
    """
    job = get_batch().new_job(f'Joint Segmentation', job_attrs | {'tool': 'gatk'})
    job.declare_resource_group(
        output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    job.image(image_path('gatk_gcnv'))
    job.memory('8Gi')
    job.command(
        f"""
    set -e
    gatk --java-options "-Xmx6000m" JointGermlineCNVSegmentation \\
    -R {reference.base} -O {job.output["vcf.gz"]} {' -V '.join(segment_vcfs)} --model-call-intervals {intervals} -ped {pedigree}
    tabix {job.output["vcf.gz"]}
    """
    )
    return job, job.output


def run_joint_segmentation(
    segment_vcfs: list[ResourceFile],
    pedigree: str,
    intervals: str,
    tmp_prefix: str,
    output_path: Path,
    job_attrs: dict[str, str] | None = None,
) -> list[Job]:
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates

    Returns:

    """

    if can_reuse(output_path):
        return []

    jobs = []

    pedigree_in_batch = get_batch().read_input(pedigree)
    intervals_in_batch = get_batch().read_input(intervals)

    # find the number of samples to shove into each scatter block
    sams_per_block = get_config()['workflow']['num_samples_per_scatter_block']

    reference = fasta_res_group(get_batch())

    chunked_vcfs = []

    # if we have more samples to process than the block size, condense
    if len(segment_vcfs) > sams_per_block:
        for subchunk_index, chunk_vcfs in enumerate(
            chunks(segment_vcfs, sams_per_block)
        ):
            # create a new job for each chunk
            # read these files into this batch
            print(chunk_vcfs)
            local_vcfs = [
                get_batch().read_input_group(vcf=vcf, index=f'{vcf}.tbi')['vcf']
                for vcf in chunk_vcfs
            ]
            job, vcf_group = joint_segment_vcfs(
                local_vcfs,
                pedigree=pedigree_in_batch,
                reference=reference,
                intervals=intervals_in_batch,
                job_attrs=job_attrs or {} | {'title': f'sub-chunk_{subchunk_index}'},
            )
            chunked_vcfs.append(vcf_group)
            get_batch().write_output(
                vcf_group, f'{tmp_prefix}/subchunk_{subchunk_index}'
            )
            jobs.append(job)

    # else, just read those into the batch
    else:
        chunked_vcfs = [
            get_batch().read_input_group(vcf=vcf, index=f'{vcf}.tbi').vcf
            for vcf in segment_vcfs
        ]

    # second round of condensing output
    job, vcf_group = joint_segment_vcfs(
        chunked_vcfs,
        pedigree=pedigree_in_batch,
        reference=reference,
        intervals=intervals_in_batch,
        job_attrs=job_attrs or {} | {'title': f'all-chunks'},
    )
    jobs.append(job)

    # write the final output file
    get_batch().write_output(vcf_group, f'{str(output_path).removesuffix(".vcf.gz")}')
    return jobs


def merge_calls(
    sg_vcfs: list[str], docker_image: str, job_attrs: dict[str, str], output_path: Path
):
    """
    This job will run a fast simple merge on per-SGID call files
    It then throws in a python script to add in two additional header lines
    and edit the SVLEN and SVTYPE attributes into each row

    Args:
        sg_vcfs (list[str]): paths to all individual VCFs
        docker_image (str): docker image to use
        job_attrs (dict): any params to atach to the job
        output_path (Path): path to the final merged VCF
    """

    if can_reuse(output_path):
        return None

    assert sg_vcfs, 'No VCFs to merge'

    merge_job = get_batch().new_job(
        'Merge gCNV calls', job_attrs | {'tool': 'bcftools'}
    )
    merge_job.image(docker_image)

    # this should be made reactive, in case we scale past 10GB
    merge_job.storage('10Gi')

    batch_vcfs = []
    for each_vcf in sg_vcfs:
        batch_vcfs.append(
            get_batch().read_input_group(
                **{
                    'vcf.gz': each_vcf,
                    'vcf.gz.tbi': f'{each_vcf}.tbi',
                }
            )['vcf.gz']
        )

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: compression level
    merge_job.command(
        f'bcftools merge {" ".join(batch_vcfs)} -Oz -o {merge_job.tmp_vcf} --threads 4 -m all -0'
    )

    # create a python job to do the file content updates
    pyjob = get_batch().new_python_job('Update VCF content')
    pyjob.storage('10Gi')
    pyjob.call(update_vcf_attributes, merge_job.tmp_vcf, pyjob.output)

    # a third job just to tidy up
    third_job = get_batch().new_job('bgzip and tabix')
    third_job.image(docker_image)
    third_job.declare_resource_group(
        output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'}
    )
    third_job.command(f'bgzip -c {pyjob.output} > {third_job.output["vcf.bgz"]}')
    third_job.command(f'tabix {third_job.output["vcf.bgz"]}')

    # dependency setting between jobs should be implicit due to temp file passing

    # get the output root to write to
    output_no_suffix = str(output_path).removesuffix('.vcf.bgz')
    get_batch().write_output(third_job.output, output_no_suffix)
    return [merge_job, pyjob, third_job]


def update_vcf_attributes(input_tmp: str, output_file: str):
    """
    A Python method to call as a PythonJob, edits content of the VCF
    - Add 2 new INFO sections in the header, SVTYPE and SVLEN
    - Use the Alt-allele post splitting to find DUP/DEL for each line
    - Use the END value (in INFO) to determine CNV Length (SVLEN)
    - Update the ID field to be unique for each line
    - Expand the INFO in each line
    - write the file back out to the specified path

    Args:
        input_tmp (str): path to temp file generated by merging
        output_file (str): path to write uncompressed edited version to
    """
    import gzip

    headers = []
    others = []

    # read the merged gVCF
    with gzip.open(input_tmp, 'rt') as f:
        for line in f:
            # don't alter current header lines
            if line.startswith('#'):
                headers.append(line)
                # but do insert additional INFO field lines
                if line.startswith('##INFO=<ID=END'):
                    headers.extend(
                        [
                            '##INFO=<ID=SVTYPE,Number=.,Type=String,Description="SV Type">\n',
                            '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV Length">\n',
                        ]
                    )
            # for non-header lines
            else:
                # split on tabs
                l_split = line.split('\t')
                original_start = int(l_split[1])

                # e.g. END=12345
                # this will be added back in upon export
                original_end = l_split[7]

                # steal the END integer (only current INFO field, int)
                end_int = int(original_end.removeprefix('END='))

                # e.g. <DEL> -> DEL
                alt_allele = l_split[4][1:-1]

                # grab the original ID
                original_id = l_split[2]

                # make this unique after splitting (include alt allele)
                l_split[2] = f'{original_id}_{alt_allele}'

                # update the INFO field with Length and Type (DUP/DEL, not "CNV")
                l_split[
                    7
                ] = f'SVTYPE={alt_allele};SVLEN={end_int - original_start};{original_end}'

                # put it together and what have you got?
                # bippidy boppidy boo
                others.append('\t'.join(l_split))

    with open(output_file, 'w') as f:
        f.writelines(headers)
        f.writelines(others)
