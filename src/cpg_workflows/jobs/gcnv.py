"""
Jobs that implement GATK-gCNV.
"""

from collections.abc import Iterable

from hailtop.batch.job import BashJob, Job
from hailtop.batch.resource import JobResourceFile, Resource, ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import config_retrieve, get_config, image_path
from cpg_utils.hail_batch import command, fasta_res_group, get_batch, query_command
from cpg_workflows.filetypes import CramPath
from cpg_workflows.query_modules import seqr_loader, seqr_loader_cnv
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.scripts import upgrade_ped_with_inferred
from cpg_workflows.utils import can_reuse, chunks


def upgrade_ped_file(local_ped: ResourceFile, new_output: str, aneuploidies: str, ploidy_tar: str):
    """
    Update the default Pedigree with the inferred ploidy information
    update a ped file, and

    Args:
        local_ped ():
        new_output ():
        aneuploidies (str): where to write identified aneuploidies
        ploidy_tar ():
    """

    j = get_batch().new_bash_job('Upgrade PED file with inferred Ploidy')
    j.image(config_retrieve(['workflow', 'driver_image']))

    # path to the python script
    script_path = upgrade_ped_with_inferred.__file__.removeprefix('/production-pipelines')
    j.command(f'tar -xf {ploidy_tar} -C .')  # creates the folder ploidy-calls
    j.command(f'python3 {script_path} {local_ped} {j.output} {j.aneuploidies} ploidy-calls')
    get_batch().write_output(j.output, new_output)
    get_batch().write_output(j.aneuploidies, aneuploidies)
    return j


def prepare_intervals(job_attrs: dict[str, str], output_paths: dict[str, Path]) -> Job:
    j = get_batch().new_job(
        'Prepare intervals',
        job_attrs | {'tool': 'gatk PreprocessIntervals/AnnotateIntervals'},
    )
    j.image(image_path('gatk_gcnv'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(get_batch())

    exclude_intervals = get_config()['workflow'].get('exclude_intervals', [])
    exclude_intervals_args = ' '.join([f'--exclude-intervals {i}' for i in exclude_intervals])

    if sequencing_type == 'exome':
        intervals = get_batch().read_input(get_config()['workflow'].get('intervals_path'))
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
    j = get_batch().new_job('Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'})
    j.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(get_batch())

    j.declare_resource_group(
        counts={
            'counts.tsv.gz': '{root}.counts.tsv.gz',
            'counts.tsv.gz.tbi': '{root}.counts.tsv.gz.tbi',
        },
    )

    assert isinstance(j.counts, ResourceGroup)

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(j)

    cmd = f"""
    gatk --java-options "{job_res.java_mem_options()}" CollectReadCounts \\
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
            },
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

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(j)

    counts_input_args = _counts_input_args(counts_paths)
    cmd = ''

    if can_reuse(output_paths['filtered']):
        # Remove 'filtered' entry from output_paths so we don't write_output() it later
        filtered: ResourceFile = get_batch().read_input(str(output_paths.pop('filtered')))
    else:
        preprocessed_intervals = get_batch().read_input(str(preprocessed_intervals_path))
        annotated_intervals = get_batch().read_input(str(annotated_intervals_path))

        cmd += f"""
        gatk --java-options "{job_res.java_mem_options()}" FilterIntervals \\
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
    gatk --java-options "{job_res.java_mem_options()}" DetermineGermlineContigPloidy \\
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

        # set highmem resources for this job
        job_res = HIGHMEM.request_resources(ncpu=8, mem_gb=52, storage_gb=10)
        job_res.set_to_job(j)

        cmd = f"""
        tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR

        {select_cmd} < {filtered_intervals} > {j.shard_intervals}

        gatk --java-options "{job_res.java_mem_options()}" GermlineCNVCaller \\
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


def postprocess_calls(
    ploidy_calls_path: Path,
    shard_paths: dict[str, Path],
    sample_index: int,
    job_attrs: dict[str, str],
    output_prefix: str,
    clustered_vcf: str | None = None,
    intervals_vcf: str | None = None,
    qc_file: str | None = None,
) -> Job:
    if any([clustered_vcf, intervals_vcf, qc_file]):
        assert all([clustered_vcf, intervals_vcf, qc_file]), [clustered_vcf, intervals_vcf, qc_file]

    j = get_batch().new_job('Postprocess gCNV calls', job_attrs | {'tool': 'gatk PostprocessGermlineCNVCalls'})
    j.image(image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=20)
    job_res.set_to_job(j)

    reference = fasta_res_group(get_batch())

    ploidy_calls_tarball = get_batch().read_input(str(ploidy_calls_path))

    gcp_related_commands = [f'tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR/inputs']

    model_shard_args = ''
    calls_shard_args = ''

    # import with hail batch instead, then unpack in the container. No need to depend on a gcloud/gsutil install
    # https://batch.hail.populationgenomics.org.au/batches/454143/jobs/1
    for name, path in [(shard, shard_paths[shard]) for shard in shard_basenames()]:
        shard_tar = get_batch().read_input(str(path))
        gcp_related_commands.append(f'tar -xzf {shard_tar} -C $BATCH_TMPDIR/inputs')
        model_shard_args += f' --model-shard-path $BATCH_TMPDIR/inputs/{name}-model'
        calls_shard_args += f' --calls-shard-path $BATCH_TMPDIR/inputs/{name}-calls'

    j.command(command(gcp_related_commands, setup_gcp=True))

    allosomal_contigs = get_config()['workflow'].get('allosomal_contigs', [])
    allosomal_contigs_args = ' '.join([f'--allosomal-contig {c}' for c in allosomal_contigs])

    # declare all output files in advance
    j.declare_resource_group(
        output={
            'intervals.vcf.gz': '{root}/intervals.vcf.gz',
            'intervals.vcf.gz.tbi': '{root}/intervals.vcf.gz.tbi',
            'segments.vcf.gz': '{root}/segments.vcf.gz',
            'segments.vcf.gz.tbi': '{root}/segments.vcf.gz.tbi',
            'ratios.tsv': '{root}/ratios.tsv',
        },
    )

    extra_args = ''
    if clustered_vcf:
        assert isinstance(intervals_vcf, str)
        local_clusters = get_batch().read_input_group(vcf=clustered_vcf, index=f'{clustered_vcf}.tbi').vcf
        local_intervals = get_batch().read_input_group(vcf=intervals_vcf, index=f'{intervals_vcf}.tbi').vcf
        extra_args += f"""--clustered-breakpoints {local_clusters} \\
         --input-intervals-vcf {local_intervals} \\
          -R {reference.base}
        """

    j.command(
        f"""
    OUTS=$(dirname {j.output['intervals.vcf.gz']})
    BATCH_OUTS=$(dirname $OUTS)
    mkdir $OUTS
    gatk --java-options "{job_res.java_mem_options()}" PostprocessGermlineCNVCalls \\
      --sequence-dictionary {reference.dict} {allosomal_contigs_args} \\
      --contig-ploidy-calls $BATCH_TMPDIR/inputs/ploidy-calls \\
      {model_shard_args} {calls_shard_args} \\
      --sample-index {sample_index} \\
      --output-genotyped-intervals {j.output['intervals.vcf.gz']} \\
      --output-genotyped-segments {j.output['segments.vcf.gz']} \\
      --output-denoised-copy-ratios {j.output['ratios.tsv']} {extra_args}
    """,
    )

    # index the output VCFs
    j.command(f'tabix -f {j.output["intervals.vcf.gz"]}')  # type: ignore
    j.command(f'tabix -f {j.output["segments.vcf.gz"]}')  # type: ignore

    if clustered_vcf:
        assert isinstance(qc_file, str)
        max_events = get_config()['workflow']['gncv_max_events']
        max_pass_events = get_config()['workflow']['gncv_max_pass_events']

        # do some additional QC to determine pass/fail
        j.command(
            f"""
        #use awk instead of grep - grep returning no lines causes a pipefailure
        NUM_SEGMENTS=$(zcat {j.output['segments.vcf.gz']} | awk '!/^#/ && !/0\/0/ && !/\t0:1:/ {{count++}} END {{print count}}')
        NUM_PASS_SEGMENTS=$(zcat {j.output['segments.vcf.gz']} | awk '!/^#/ && !/0\/0/ && !/\t0:1:/ && /PASS/ {{count++}} END {{print count}}')
        if [ $NUM_SEGMENTS -lt {max_events} ]; then
            if [ $NUM_PASS_SEGMENTS -lt {max_pass_events} ]; then
              echo "PASS" >> {j.qc_file}
            else
              echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" >> {j.qc_file}
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> {j.qc_file}
        fi
        cat {j.qc_file}
        """,
        )
        get_batch().write_output(j.qc_file, qc_file)

    get_batch().write_output(j.output, output_prefix)

    return j


def joint_segment_vcfs(
    segment_vcfs: list[ResourceFile],
    pedigree: ResourceFile,
    reference: ResourceGroup,
    intervals: ResourceFile,
    title: str,
    job_attrs: dict,
) -> tuple[BashJob, Resource]:
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates

    Returns:
        the job that does the work, and the resulting resource group of VCF & index
    """
    job = get_batch().new_bash_job(f'Joint Segmentation {title}', job_attrs | {'tool': 'gatk'})
    job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    job.image(image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(job)

    vcf_string = ''
    for each_vcf in segment_vcfs:
        vcf_string += f' -V {each_vcf}'

    # this already creates a tabix index
    job.command(
        f"""
    set -e
    gatk --java-options "{job_res.java_mem_options()}" JointGermlineCNVSegmentation \\
        -R {reference.base} \\
        -O {job.output["vcf.gz"]} \\
        {vcf_string} \\
        --model-call-intervals {intervals} \\
        -ped {pedigree}
    """,
    )
    return job, job.output


def run_joint_segmentation(
    segment_vcfs: list[str],
    pedigree: str,
    intervals: str,
    tmp_prefix: str,
    output_path: Path,
    job_attrs: dict[str, str] | None = None,
) -> list[BashJob]:
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates
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
    # this is done by calling an intermediate round of VCF segmenting
    if len(segment_vcfs) > sams_per_block:
        for subchunk_index, chunk_vcfs in enumerate(chunks(segment_vcfs, sams_per_block)):
            # create a new job for each chunk
            # read these files into this batch
            local_vcfs = [get_batch().read_input_group(vcf=vcf, index=f'{vcf}.tbi')['vcf'] for vcf in chunk_vcfs]
            job, vcf_group = joint_segment_vcfs(
                local_vcfs,
                pedigree=pedigree_in_batch,
                reference=reference,
                intervals=intervals_in_batch,
                job_attrs=job_attrs or {} | {'title': f'sub-chunk_{subchunk_index}'},
                title=f'sub-chunk_{subchunk_index}',
            )
            chunked_vcfs.append(vcf_group['vcf.gz'])  # type: ignore
            get_batch().write_output(vcf_group, f'{tmp_prefix}/subchunk_{subchunk_index}')
            jobs.append(job)

    # else, all vcf files into the batch
    else:
        chunked_vcfs = [get_batch().read_input_group(vcf=vcf, index=f'{vcf}.tbi').vcf for vcf in segment_vcfs]

    # second round of condensing output - produces one single file
    job, vcf_group = joint_segment_vcfs(
        chunked_vcfs,
        pedigree=pedigree_in_batch,
        reference=reference,
        intervals=intervals_in_batch,
        job_attrs=job_attrs or {} | {'title': 'all-chunks'},
        title='all-chunks',
    )
    jobs.append(job)

    # write the final output file group (VCF & index)
    get_batch().write_output(vcf_group, f'{str(output_path).removesuffix(".vcf.gz")}')
    return jobs


def trim_sex_chromosomes(sgid: str, sg_vcf: str, no_xy_vcf: str, job_attrs: dict[str, str]) -> BashJob:
    """
    Create a BCFtools job to trim chrX and chrY from this VCF

    Args:
        sgid (str): the SG ID, used in naming the job
        sg_vcf (str): the path to the input VCF
        no_xy_vcf (str): the path to the output VCF
        job_attrs ():

    Returns:
        the job which generates the output file
    """
    job = get_batch().new_bash_job(f'remove sex chromosomes from {sgid}', job_attrs | {'tool': 'bcftools'})
    job.image(image_path('bcftools'))
    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
    localised_vcf = get_batch().read_input_group(**{'vcf.gz': sg_vcf, 'vcf.gz.tbi': f'{sg_vcf}.tbi'})['vcf.gz']
    autosomes = ' '.join([f'chr{i}' for i in range(1, 23)])
    job.command('set -euo pipefail')
    job.command(f'bcftools view {localised_vcf} {autosomes} | bgzip -c > {job.output["vcf.bgz"]}')  # type: ignore
    job.command(f'tabix {job.output["vcf.bgz"]}')  # type: ignore
    get_batch().write_output(job.output, no_xy_vcf.removesuffix('.vcf.bgz'))
    return job


def merge_calls(sg_vcfs: list[str], docker_image: str, job_attrs: dict[str, str], output_path: Path):
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

    merge_job = get_batch().new_job('Merge gCNV calls', job_attrs | {'tool': 'bcftools'})
    merge_job.image(docker_image)

    # this should be made reactive, in case we scale past 10GB
    merge_job.storage('10Gi')

    batch_vcfs = []
    for each_vcf in sg_vcfs:
        batch_vcfs.append(
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz'],
        )

    # option breakdown:
    # -Oz: bgzip output
    # -o: output file
    # --threads: number of threads to use
    # -m: merge strategy
    # -0: compression level
    merge_job.command(f'bcftools merge {" ".join(batch_vcfs)} -Oz -o {merge_job.tmp_vcf} --threads 4 -m all -0')

    # now normlise the result, splitting multiallelics
    merge_job.command(f'bcftools norm -m -any {merge_job.tmp_vcf} | bgzip -c > {merge_job.tmp_vcf_split}')

    # create a python job to do the file content updates
    pyjob = get_batch().new_python_job('Update VCF content')
    pyjob.storage('10Gi')
    pyjob.call(update_vcf_attributes, merge_job.tmp_vcf_split, pyjob.output)

    # a third job just to tidy up
    third_job = get_batch().new_job('bgzip and tabix')
    third_job.image(docker_image)
    third_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
    third_job.command(f'bcftools view {pyjob.output} | bgzip -c > {third_job.output["vcf.bgz"]}')  # type: ignore
    third_job.command(f'tabix {third_job.output["vcf.bgz"]}')  # type: ignore

    # dependency setting between jobs should be implicit due to temp file passing

    # get the output root to write to
    output_no_suffix = str(output_path).removesuffix('.vcf.bgz')
    get_batch().write_output(third_job.output, output_no_suffix)
    return [merge_job, pyjob, third_job]


def update_vcf_attributes(input_tmp: str, output_file: str):
    """
    A Python method to call as a PythonJob, edits content of the VCF
    - INFO.SVLEN exists in header, but not populated on variant rows
    - Use the END value (in INFO) to determine CNV Length (SVLEN)
    - Update the ID field to be appropriate for each line
    - Expand the INFO in each line
    - write the file back out to the specified path

    Read and write are both done compressed

    Args:
        input_tmp (str): path to temp file generated by merging
        output_file (str): path to write uncompressed edited version to
    """
    import gzip

    # read the merged gVCF
    with gzip.open(input_tmp, 'rt') as f:
        with gzip.open(output_file, 'wt') as f_out:
            for line in f:
                # don't alter current header lines
                if line.startswith('#'):
                    f_out.write(line)
                    continue

                # for non-header lines, split on tabs
                l_split = line.split('\t')

                # don't bother with null/WT/missing alleles
                if l_split[4] == '.':
                    continue

                original_start = int(l_split[1])

                # e.g. AN_Orig=61;END=56855888;SVTYPE=DUP
                original_info: dict[str, str] = dict(el.split('=') for el in l_split[7].split(';'))

                # steal the END integer
                end_int = int(original_info['END'])

                # e.g. <DEL> -> DEL
                alt_allele = l_split[4][1:-1]

                # replace compound ID original ID
                chrom = l_split[0]
                start = l_split[1]

                # make this unique after splitting (include alt allele)
                l_split[2] = f'CNV_{chrom}_{start}_{end_int}_{alt_allele}'

                # update the INFO field with Length and Type (DUP/DEL, not "CNV")
                original_info['SVLEN'] = str(end_int - original_start)
                l_split[7] = ';'.join(f'{k}={v}' for k, v in original_info.items())

                # put it together and what have you got?
                # bippidy boppidy boo
                f_out.write('\t'.join(l_split))


def annotate_dataset_jobs_cnv(
    mt_path: Path,
    sgids: list[str],
    out_mt_path: Path,
    tmp_prefix: Path,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Split mt by dataset and annotate dataset-specific fields (only for those datasets
    that will be loaded into Seqr).
    """
    assert sgids
    sgids_list_path = tmp_prefix / 'sgid-list.txt'
    if not get_config()['workflow'].get('dry_run', False):
        with sgids_list_path.open('w') as f:
            f.write(','.join(sgids))

    subset_mt_path = tmp_prefix / 'cohort-subset.mt'
    subset_j: Job | None = None
    if not subset_mt_path.exists():
        subset_j = get_batch().new_job('subset cohort to dataset', (job_attrs or {}) | {'tool': 'hail query'})
        subset_j.image(image_path('cpg_workflows'))
        subset_j.command(
            query_command(
                seqr_loader,
                seqr_loader.subset_mt_to_samples.__name__,
                str(mt_path),
                sgids,
                str(subset_mt_path),
                setup_gcp=True,
            ),
        )

    annotate_j = get_batch().new_job('annotate dataset', (job_attrs or {}) | {'tool': 'hail query'})
    annotate_j.image(image_path('cpg_workflows'))
    annotate_j.command(
        query_command(
            seqr_loader_cnv,
            seqr_loader_cnv.annotate_dataset_gcnv.__name__,
            str(subset_mt_path),
            str(out_mt_path),
            setup_gcp=True,
        ),
    )
    if subset_j:
        annotate_j.depends_on(subset_j)
        return [subset_j, annotate_j]
    else:
        return [annotate_j]
