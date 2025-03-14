"""
Create and run jobs relating to VQSR for the RD combiner
"""

from functools import lru_cache

from hailtop.batch.job import Job
from hailtop.batch.resource import Resource, ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.vcf import quick_and_easy_bcftools_concat
from cpg_workflows.jobs.vqsr import (
    INDEL_ALLELE_SPECIFIC_FEATURES,
    INDEL_RECALIBRATION_TRANCHE_VALUES,
    SNP_ALLELE_SPECIFIC_FEATURES,
    SNP_RECALIBRATION_TRANCHE_VALUES,
)
from cpg_workflows.resources import HIGHMEM, STANDARD
from cpg_workflows.utils import (
    VCF_GZ,
    VCF_GZ_TBI,
    can_reuse,
    chunks,
    generator_chunks,
    get_all_fragments_from_manifest,
    get_logger,
)

TRAINING_PER_JOB: int = config_retrieve(['rd_combiner', 'vqsr_training_fragments_per_job'], 100)
RECALIBRATION_PER_JOB: int = config_retrieve(['rd_combiner', 'vqsr_apply_fragments_per_job'], 60)
INDEL_RECAL_DISC_SIZE: int = config_retrieve(['rd_combiner', 'indel_recal_disc_size'], 20)
SNPS_RECAL_DISC_SIZE: int = config_retrieve(['rd_combiner', 'snps_recal_disc_size'], 20)
SNPS_GATHER_DISC_SIZE: int = config_retrieve(['rd_combiner', 'snps_gather_disc_size'], 10)


@lru_cache(1)
def get_localised_resources_for_vqsr() -> dict[str, ResourceGroup]:
    """
    get the resources required for VQSR, once per run
    Returns:
        the dictionary of resources and their names
    """

    return {
        key: get_batch().read_input_group(
            base=reference_path(f'broad/{key}_vcf'),
            index=reference_path(f'broad/{key}_vcf_index'),
        )
        for key in [
            'axiom_poly',
            'dbsnp',
            'hapmap',
            'mills',
            'omni',
            'one_thousand_genomes',
        ]
    }


def train_vqsr_indels(sites_only_vcf: str, output_prefix: str, job_attrs: dict) -> Job:
    """
    Train VQSR indels on the sites-only VCF
    Args:
        sites_only_vcf (str): path to the whole-genome sites-only VCF
        output_prefix (str): root path to write the output files to
        job_attrs (dict): job attributes
    Returns:
        a Job object
    """
    resources = get_localised_resources_for_vqsr()
    siteonly_vcf = get_batch().read_input_group(
        **{
            VCF_GZ: str(sites_only_vcf),
            VCF_GZ_TBI: str(sites_only_vcf) + '.tbi',
        },
    )
    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs.
    """
    indel_recalibrator_j = get_batch().new_job(
        'TrainVqsrIndelModelOnCombinerData',
        job_attrs | {'tool': 'gatk VariantRecalibrator'},
    )
    indel_recalibrator_j.image(image_path('gatk'))
    indel_recalibrator_j.command('set -euo pipefail')

    # We run it for the entire dataset in one job, so can take an entire instance.
    instance_fraction = 1
    res = HIGHMEM.set_resources(indel_recalibrator_j, fraction=instance_fraction, storage_gb=INDEL_RECAL_DISC_SIZE)

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in INDEL_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [f'-an {v}' for v in INDEL_ALLELE_SPECIFIC_FEATURES],
    )

    # delclare a resource group for the output
    indel_recalibrator_j.declare_resource_group(
        output={
            'recal': '{root}.recal',
            'recal.idx': '{root}.recal.idx',
            'tranches': '{root}.tranches',
        },
    )
    indel_recalibrator_j.command(
        f"""
        gatk --java-options \
          "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
          VariantRecalibrator \\
          -V {siteonly_vcf['vcf.gz']} \\
          -O {indel_recalibrator_j.output.recal} \\
          --tranches-file {indel_recalibrator_j.output.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          --use-allele-specific-annotations \\
          --max-gaussians 4 \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {resources['mills'].base} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {resources['axiom_poly'].base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {resources['dbsnp'].base}
        """,
    )
    get_batch().write_output(indel_recalibrator_j.output, output_prefix)
    return indel_recalibrator_j


def train_vqsr_snps(sites_only_vcf: str, snp_model: str, job_attrs: dict):
    """
    Train VQSR indels on the sites-only VCF
    Args:
        sites_only_vcf ():
        snp_model ():
        job_attrs ():
    """
    resources = get_localised_resources_for_vqsr()
    siteonly_vcf = get_batch().read_input_group(
        **{
            VCF_GZ: sites_only_vcf,
            VCF_GZ_TBI: sites_only_vcf + '.tbi',
        },
    )

    snp_recalibrator_j = get_batch().new_job(
        'TrainVqsrSnpModelOnCombinerData',
        job_attrs | {'tool': 'gatk VariantRecalibrator'},
    )
    snp_recalibrator_j.image(image_path('gatk'))

    # We run it for the entire dataset in one job, so can take an entire instance.
    res = HIGHMEM.set_resources(snp_recalibrator_j, fraction=1, storage_gb=SNPS_RECAL_DISC_SIZE)

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [f'-an {v}' for v in SNP_ALLELE_SPECIFIC_FEATURES],
    )
    snp_recalibrator_j.command(
        f"""set -euo pipefail
        gatk --java-options \
          "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
          VariantRecalibrator \\
          -V {siteonly_vcf['vcf.gz']} \\
          -O {snp_recalibrator_j.recalibration} \\
          --tranches-file {snp_recalibrator_j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          --use-allele-specific-annotations \\
          --sample-every-Nth-variant 10 \\
          --output-model {snp_recalibrator_j.model_file} \\
          --max-gaussians 6 \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {resources['hapmap'].base} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {resources['omni'].base} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {resources['one_thousand_genomes'].base} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {resources['dbsnp'].base}
          """,
    )

    get_batch().write_output(snp_recalibrator_j.model_file, snp_model)
    return snp_recalibrator_j


def train_vqsr_snp_tranches(
    manifest_file: Path,
    snp_model_path: str,
    output_path: str,
    temp_path: Path,
    job_attrs: dict,
) -> list[Job]:
    """
    train vqsr tranches for SNPs, in a scattered manner
    Args:
        manifest_file (Path): Path object to the manifest file
        snp_model_path (str): path to the SNP model trained in the previous stage
        output_path (str): path to write the tranches to
        temp_path (Path): path to the temp directory
        job_attrs (dict): job attributes

    Returns:
        all jobs required to get where we're going
    """

    vcf_resources = get_all_fragments_from_manifest(manifest_file)

    fragment_count = len(vcf_resources)
    snps_recal_paths = [temp_path / f'snp_{i}.recal' for i in range(fragment_count)]
    snps_tranches_paths = [temp_path / f'snp_{i}.tranches' for i in range(fragment_count)]

    snp_model_in_batch = get_batch().read_input(snp_model_path)

    resources = get_localised_resources_for_vqsr()

    # the list of all jobs (per-fragment, and the summary)
    scatter_jobs: list[Job] = []

    # to hold the resulting parts, as Resources
    snp_tranche_fragments: list[ResourceFile] = []

    # if we start this as -1, we can increment at the start of the loop, making the index counting easier to track
    vcf_counter = -1

    # iterate over all fragments, but in chunks of FRAGMENTS_PER_JOB
    for chunk_counter, fragment_chunk in enumerate(chunks(vcf_resources, TRAINING_PER_JOB)):
        # NB this VQSR training stage is scattered, and we have a high number of very small VCF fragments
        # 99.5% of the time and cost of this task was pulling the docker image and loading reference data
        # the actual work took 5 seconds at negligible cost. Instead of running one job per VCF fragment,
        # we can batch fragments into fewer jobs, each stacking multiple fragments together but only requiring
        # one block of reference data to be loaded.
        # candidate splitting is 100-fragments-per-job, for a ~99% cost saving

        chunk_job = get_batch().new_job(f'TrainVqsrSnpTranches, Chunk {chunk_counter}', job_attrs)
        chunk_job.image(image_path('gatk'))
        chunk_job.command('set -euo pipefail')

        # add this job to the list of scatter jobs
        scatter_jobs.append(chunk_job)

        res = STANDARD.set_resources(chunk_job, ncpu=4, storage_gb=SNPS_GATHER_DISC_SIZE)

        # iterate over the fragment VCF resource groups
        for vcf_resource in fragment_chunk:
            vcf_counter += 1

            snps_recal_path = snps_recal_paths[vcf_counter]
            snps_tranche_path = snps_tranches_paths[vcf_counter]

            if can_reuse(snps_recal_path) and can_reuse(snps_tranche_path):
                get_logger().info(f'Reusing {snps_recal_path} and {snps_tranche_path}')
                snp_tranche_fragments.append(get_batch().read_input(str(snps_tranche_path)))
                continue

            tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
            an_cmdl = ' '.join(
                [f'-an {v}' for v in SNP_ALLELE_SPECIFIC_FEATURES],
            )

            # create a counter string to uniquely reference all job outputs
            counter_string = str(vcf_counter)

            # create a resource group for the recalibration output and its index
            chunk_job.declare_resource_group(
                **{
                    counter_string: {
                        'recal': '{root}.recal',
                        'recal.idx': '{root}.recal.idx',
                        'tranches': '{root}.tranches',
                    },
                },
            )

            # the mv command here is because the input VCF is a .bgz instead of .vcf.bgz
            # so GATK can't tell what type of file it is
            chunk_job.command(
                f"""
                MODEL_REPORT={snp_model_in_batch}
                mv {vcf_resource[VCF_GZ]} input.vcf.bgz
                mv {vcf_resource[VCF_GZ_TBI]} input.vcf.bgz.tbi
                gatk --java-options \
                  "{res.java_mem_options()} {res.java_gc_thread_options()}" \\
                  VariantRecalibrator \\
                  --verbosity WARNING \\
                  --QUIET \\
                  -V input.vcf.bgz \\
                  -O {chunk_job[counter_string]['recal']} \\
                  --tranches-file {chunk_job[counter_string]['tranches']} \\
                  --trust-all-polymorphic \\
                  {tranche_cmdl} \\
                  {an_cmdl} \\
                  -mode SNP \\
                  --use-allele-specific-annotations \\
                  --input-model {snp_model_in_batch} \\
                  --output-tranches-for-scatter \\
                  --max-gaussians 6 \\
                  -resource:hapmap,known=false,training=true,truth=true,prior=15 {resources['hapmap'].base} \\
                  -resource:omni,known=false,training=true,truth=true,prior=12 {resources['omni'].base} \\
                  -resource:1000G,known=false,training=true,truth=false,prior=10 {resources['one_thousand_genomes'].base} \\
                  -resource:dbsnp,known=true,training=false,truth=false,prior=7 {resources['dbsnp'].base}
                touch {chunk_job[counter_string]['recal.idx']}
                """,
            )

            # write the results out to GCP
            get_batch().write_output(chunk_job[counter_string], str(temp_path / f'snp_{vcf_counter}'))
            snp_tranche_fragments.append(chunk_job[counter_string].tranches)

    # one final job to write the success indicator
    final_job = get_batch().new_bash_job('Completion message', job_attrs)
    final_job.image(config_retrieve(['workflow', 'driver_image']))
    final_job.command(f'echo "All tranches trained" > {final_job.output}')
    final_job.depends_on(*scatter_jobs)
    scatter_jobs.append(final_job)
    get_batch().write_output(final_job.output, output_path)
    return scatter_jobs


def gather_tranches(manifest_file: Path, temp_path: Path, output_path: str, job_attrs: dict) -> Job:
    """
    The previous approach ran into hard limits on the size of the batch spec
    There was too much metadata around which resource groups the tranches were part of etc etc
    Splitting out data generating from data aggregating should hopefully help

    Args:
        manifest_file (Path): path to the manifest file
        temp_path (Path): path to the temp directory (same as previous stage)
        output_path (str): path to write the tranches aggregate to
        job_attrs (dict): job attributes
    """

    vcf_resources = get_all_fragments_from_manifest(manifest_file)
    snp_tranche_paths = [
        get_batch().read_input(str(temp_path / f'snp_{i}.tranches')) for i in range(len(vcf_resources))
    ]

    gather_tranches_j = get_batch().new_job('GatherTrainedVqsrSnpTranches', job_attrs | {'tool': 'gatk GatherTranches'})
    gather_tranches_j.image(image_path('gatk'))
    res = STANDARD.set_resources(gather_tranches_j, ncpu=2, storage_gb=SNPS_GATHER_DISC_SIZE)

    inputs_cmdl = ' '.join([f'--input {t}' for t in snp_tranche_paths])
    gather_tranches_j.command(
        f"""set -euo pipefail
    gatk --java-options "{res.java_mem_options()}" \\
      GatherTranches \\
      --mode SNP \\
      {inputs_cmdl} \\
      --output {gather_tranches_j.out_tranches}""",
    )

    get_batch().write_output(gather_tranches_j.out_tranches, output_path)
    return gather_tranches_j


def apply_snp_vqsr_to_fragments(
    manifest_file: Path,
    tranche_file: str,
    temp_path: Path,
    output_path: str,
    job_attrs: dict,
):
    """
    Apply SNP VQSR to the tranches
    I'm going to retry the stacking approach again to reduce job count

    Gather results into a single file

    Args:
        manifest_file (Path): path to the manifest file, locating all VCF fragments
        tranche_file ():
        temp_path (): Path to the temp from TrainVqsrSnpTranches
        output_path (str):
        job_attrs ():
    """

    vcf_resources = get_all_fragments_from_manifest(manifest_file)
    fragment_count = len(vcf_resources)

    # read all the recal fragments into the batch as ResourceGroups
    # we're creating these paths in expectation that they were written by the tranches stage
    snps_recal_resources = [
        get_batch().read_input_group(
            recal=str(temp_path / f'snp_{i}.recal'),
            idx=str(temp_path / f'snp_{i}.recal.idx'),
        )
        for i in range(fragment_count)
    ]

    tranches_in_batch = get_batch().read_input(tranche_file)

    applied_recalibration_jobs: list[Job] = []
    recalibrated_snp_vcfs: list[Resource] = []

    vcf_counter = -1
    snp_filter_level = config_retrieve(['vqsr', 'snp_filter_level'])

    for chunk_counter, vcfs_recals in enumerate(
        generator_chunks(zip(vcf_resources, snps_recal_resources), RECALIBRATION_PER_JOB),
    ):

        chunk_job = get_batch().new_bash_job(f'RunTrainedSnpVqsrOnCombinerFragments, Chunk {chunk_counter}', job_attrs)
        chunk_job.image(image_path('gatk'))

        # stores all the annotated VCFs in this chunk
        chunk_vcfs = []

        res = STANDARD.set_resources(chunk_job, ncpu=1, storage_gb=10)

        # iterate over the zipped resource groups
        for vcf_resource, recal_resource in vcfs_recals:
            vcf_counter += 1
            # used in namespacing the outputs
            counter_string = str(vcf_counter)

            # create a resource group for the recalibration output and its index
            chunk_job.declare_resource_group(
                **{
                    counter_string: {
                        VCF_GZ: '{root}.vcf.gz',
                        VCF_GZ_TBI: '{root}.vcf.gz.tbi',
                    },
                },
            )

            chunk_job.command(
                f"""
            gatk --java-options "{res.java_mem_options()}" \\
            ApplyVQSR \\
            -O {chunk_job[counter_string]['vcf.gz']} \\
            -V {vcf_resource['vcf.gz']} \\
            --recal-file {recal_resource.recal} \\
            --tranches-file {tranches_in_batch} \\
            --truth-sensitivity-filter-level {snp_filter_level} \\
            --use-allele-specific-annotations \\
            -mode SNP

            tabix -p vcf -f {chunk_job[counter_string]['vcf.gz']}
            """,
            )
            chunk_vcfs.append(chunk_job[counter_string])

        # concatenates all VCFs in this chunk
        chunk_concat_job = quick_and_easy_bcftools_concat(
            chunk_vcfs,
            storage_gb=SNPS_GATHER_DISC_SIZE,
            job_attrs=job_attrs,
        )
        chunk_concat_job.depends_on(chunk_job)
        applied_recalibration_jobs.extend([chunk_job, chunk_concat_job])
        recalibrated_snp_vcfs.append(chunk_concat_job.output)

    # now we've got all the recalibrated VCFs, we need to gather them into a single VCF
    final_gather_job = quick_and_easy_bcftools_concat(
        recalibrated_snp_vcfs,
        storage_gb=SNPS_GATHER_DISC_SIZE,
        job_attrs=job_attrs,
    )
    final_gather_job.depends_on(*applied_recalibration_jobs)
    applied_recalibration_jobs.append(final_gather_job)

    get_batch().write_output(final_gather_job.output, output_path.removesuffix('.vcf.gz'))
    return applied_recalibration_jobs


def apply_recalibration_indels(
    snp_annotated_vcf: Path,
    indel_recalibration: Path,
    indel_tranches: Path,
    output_path: Path,
    job_attrs: dict,
):
    """
    Apply indel recalibration to the annotated SNP VCF
    """

    snp_vcf_in_batch = get_batch().read_input_group(
        vcf=str(snp_annotated_vcf),
        vcf_index=f'{str(snp_annotated_vcf)}.tbi',
    )
    indel_tranches_in_batch = get_batch().read_input(str(indel_tranches))
    indel_recalibration_in_batch = get_batch().read_input_group(
        recal=str(indel_recalibration),
        recal_idx=f'{str(indel_recalibration)}.idx',
    )

    indel_recal_job = get_batch().new_bash_job(f'RunTrainedIndelVqsrOnCombinedVcf on {snp_annotated_vcf}', job_attrs)
    indel_recal_job.image(image_path('gatk'))
    res = STANDARD.set_resources(indel_recal_job, ncpu=2, storage_gb=INDEL_RECAL_DISC_SIZE)

    indel_recal_job.declare_resource_group(output={VCF_GZ: '{root}.vcf.gz', VCF_GZ_TBI: '{root}.vcf.gz.tbi'})

    filter_level = config_retrieve(['vqsr', 'indel_filter_level'])

    indel_recal_job.command(
        f"""
    gatk --java-options "{res.java_mem_options()}" \\
    ApplyVQSR \\
    --tmp-dir $BATCH_TMPDIR \\
    -O {indel_recal_job.output[VCF_GZ]} \\
    -V {snp_vcf_in_batch.vcf} \\
    --recal-file {indel_recalibration_in_batch.recal} \\
    --tranches-file {indel_tranches_in_batch} \\
    --truth-sensitivity-filter-level {filter_level} \\
    --use-allele-specific-annotations \\
    -mode INDEL
    tabix -p vcf -f {indel_recal_job.output[VCF_GZ]}
    """,
    )
    get_batch().write_output(indel_recal_job.output, str(output_path).removesuffix('.vcf.gz'))
    return indel_recal_job
