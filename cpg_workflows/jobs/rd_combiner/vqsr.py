"""
Create and run jobs relating to VQSR for the RD combiner
"""

from functools import lru_cache

from hailtop.batch.job import Job
from hailtop.batch.resource import ResourceFile, ResourceGroup

from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.vqsr import (
    SNP_ALLELE_SPECIFIC_FEATURES,
    SNP_RECALIBRATION_TRANCHE_VALUES,
    indel_recalibrator_job,
    snps_gather_tranches_job,
    snps_recalibrator_create_model_job,
)
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.utils import VCF_GZ, VCF_GZ_TBI, can_reuse, chunks

FRAGMENTS_PER_JOB: int = 80


@lru_cache(2)
def get_all_fragments_from_manifest(manifest_file: Path) -> list[ResourceGroup]:
    """
    read the manifest file, and return all the fragment resources as an ordered list
    this is a cached method as we don't want to localise every fragment once per task

    Args:
        manifest_file ():

    Returns:
        an ordered list of all the fragment VCFs and corresponding indices
    """

    resource_objects: list[ResourceGroup] = []
    manifest_folder: Path = manifest_file.parent
    with manifest_file.open() as f:
        for line in f:
            vcf_path = manifest_folder / line.strip()
            resource_objects.append(
                get_batch().read_input_group(**{VCF_GZ: vcf_path, VCF_GZ_TBI: f'{vcf_path}.tbi'}),
            )
    return resource_objects


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


def train_vqsr_indels(sites_only_vcf: str, indel_recal: str, indel_tranches: str):
    """
    Train VQSR indels on the sites-only VCF
    Args:
        sites_only_vcf ():
        indel_recal ():
        indel_tranches ():
    """
    resources = get_localised_resources_for_vqsr()
    siteonly_vcf = get_batch().read_input_group(
        **{
            VCF_GZ: str(sites_only_vcf),
            VCF_GZ_TBI: str(sites_only_vcf) + '.tbi',
        },
    )

    indel_recalibrator_j = indel_recalibrator_job(
        b=get_batch(),
        siteonly_vcf=siteonly_vcf,
        mills_resource_vcf=resources['mills'],
        axiom_poly_resource_vcf=resources['axiom_poly'],
        dbsnp_resource_vcf=resources['dbsnp'],
        disk_size=20,
        use_as_annotations=True,
        is_small_callset=False,
    )

    get_batch().write_output(indel_recalibrator_j.recalibration, indel_recal)
    get_batch().write_output(indel_recalibrator_j.recalibration_idx, indel_recal + '.idx')
    get_batch().write_output(indel_recalibrator_j.tranches, indel_tranches)
    return indel_recalibrator_j


def train_vqsr_snps(sites_only_vcf: str, snp_model: str):
    """
    Train VQSR indels on the sites-only VCF
    Args:
        sites_only_vcf ():
        snp_model ():
    """
    resources = get_localised_resources_for_vqsr()
    siteonly_vcf = get_batch().read_input_group(
        **{
            VCF_GZ: str(sites_only_vcf),
            VCF_GZ_TBI: str(sites_only_vcf) + '.tbi',
        },
    )

    snp_recalibrator_j = snps_recalibrator_create_model_job(
        b=get_batch(),
        siteonly_vcf=siteonly_vcf,
        hapmap_resource_vcf=resources['hapmap'],
        omni_resource_vcf=resources['omni'],
        one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
        dbsnp_resource_vcf=resources['dbsnp'],
        disk_size=20,
        use_as_annotations=True,
        is_small_callset=False,
    )

    get_batch().write_output(snp_recalibrator_j.model_file, snp_model)
    return snp_recalibrator_j


def train_vqsr_snp_tranches(
    manifest_file: Path,
    snp_model_path: str,
    output_path: str,
    temp_path: Path,
    job_name: str = 'TrainVqsrSnpTranches',
) -> list[Job]:
    """
    train vqsr tranches for SNPs, in a scattered manner
    Args:
        manifest_file (Path): Path object to the manifest file
        snp_model_path ():
        output_path ():
        temp_path ():

    Returns:
        all jobs required to get where we're going
    """

    vcf_resources = get_all_fragments_from_manifest(manifest_file)

    fragment_count = len(vcf_resources)
    snps_recal_paths = [temp_path / f'snp_recalibrations_{i}' for i in range(fragment_count)]
    snps_tranches_paths = [temp_path / f'snp_tranches_{i}' for i in range(fragment_count)]

    snp_model_in_batch = get_batch().read_input(snp_model_path)

    resources = get_localised_resources_for_vqsr()

    # the list of all jobs (per-fragment, and the summary)
    scatter_jobs: list[Job] = []

    # to hold the resulting parts, as Resources
    snp_tranche_fragments: list[ResourceFile] = []

    # if we start this as -1, we can increment at the start of the loop, making the index counting easier to track
    vcf_counter = -1

    chunk_counter = 0

    # iterate over all fragments, but in chunks of FRAGMENTS_PER_JOB
    for fragment_chunk in chunks(vcf_resources, FRAGMENTS_PER_JOB):
        # NB this VQSR training stage is scattered, and we have a high number of very small VCF fragments
        # 99.5% of the time and cost of this task was pulling the docker image and loading reference data
        # the actual work took 5 seconds at negligible cost. Instead of running one job per VCF fragment,
        # we can batch fragments into fewer jobs, each stacking multiple fragments together but only requiring
        # one block of reference data to be loaded.
        # candidate splitting is 100-fragments-per-job, for a ~99% cost saving

        chunk_counter += 1

        chunk_job = get_batch().new_job(f'{job_name}, Chunk {chunk_counter}')
        chunk_job.image(image_path('gatk'))

        # add this job to the list of scatter jobs
        scatter_jobs.append(chunk_job)

        res = HIGHMEM.set_resources(chunk_job, ncpu=4, storage_gb=50)

        # iterate over the fragment VCF resource groups
        for vcf_resource in fragment_chunk:
            vcf_counter += 1

            snps_recal_path = snps_recal_paths[vcf_counter]
            snps_tranche_path = snps_tranches_paths[vcf_counter]

            if can_reuse(snps_recal_path) and can_reuse(snps_tranche_path):
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
                        'recal': '{root}',
                        'recal_idx': '{root}.idx',
                        'tranches': '{root}.tranches',
                    },
                },
            )

            chunk_job.command(
                f"""
                set -euo pipefail
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
                touch {chunk_job[counter_string]['recal_idx']}
                """,
            )

            # write the results out to GCP
            get_batch().write_output(chunk_job[counter_string].recal, str(snps_recal_paths[vcf_counter]))
            get_batch().write_output(chunk_job[counter_string].recal_idx, str(snps_recal_paths[vcf_counter]) + '.idx')
            get_batch().write_output(chunk_job[counter_string].tranches, str(snps_tranches_paths[vcf_counter]))
            snp_tranche_fragments.append(chunk_job[counter_string].tranches)

    # one final job to write the success indicator
    final_job = get_batch().new_bash_job('Completion message')
    final_job.image(config_retrieve(['workflow', 'driver_image']))
    final_job.command(f'echo "All tranches trained" > {final_job.output}')
    final_job.depends_on(*scatter_jobs)
    scatter_jobs.append(final_job)
    get_batch().write_output(final_job.output, output_path)

    return scatter_jobs


def gather_tranches(manifest_file: Path, temp_path: Path, output_path: str) -> Job:
    """
    The previous approach ran into hard limits on the size of the batch spec
    There was too much metadata around which resource groups the tranches were part of etc etc
    Splitting out data generating from data aggregating should hopefully help

    Args:
        manifest_file (Path): path to the manifest file
        temp_path (Path): path to the temp directory (same as previous stage)
        output_path (str): path to write the tranches aggregate to
    """

    vcf_resources = get_all_fragments_from_manifest(manifest_file)
    snp_tranche_paths = [
        get_batch().read_input(str(temp_path / f'snp_tranches_{i}')) for i in range(len(vcf_resources))
    ]

    gather_tranches_j = snps_gather_tranches_job(
        get_batch(),
        tranches=snp_tranche_paths,
        disk_size=100,
    )
    get_batch().write_output(gather_tranches_j.out_tranches, output_path)
    return gather_tranches_j
