"""
Create and run jobs relating to VQSR for the RD combiner
"""

from functools import lru_cache

from hailtop.batch.resource import ResourceGroup
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.vqsr import indel_recalibrator_job, snps_recalibrator_create_model_job, snps_recalibrator_scattered, snps_gather_tranches_job
from cpg_workflows.utils import can_reuse


@lru_cache(100)
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
                get_batch().read_input_group(**{'vcf.bgz': vcf_path, 'vcf.bgz.tbi': f'{vcf_path}.tbi'}),
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
            'vcf.gz': str(sites_only_vcf),
            'vcf.gz.tbi': str(sites_only_vcf) + '.tbi',
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
            'vcf.gz': str(sites_only_vcf),
            'vcf.gz.tbi': str(sites_only_vcf) + '.tbi',
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
    manifest_file:Path,
    snp_model_path:str,
    output_path:str,
    temp_path:Path,
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

    # to hold the resulting parts, as ResourceFiles
    snp_tranche_fragments = []

    # iterate over all fragments, make
    for idx in range(fragment_count):
        snps_recal_path = snps_recal_paths[idx]
        snps_tranche_path = snps_tranches_paths[idx]

        if can_reuse(snps_recal_path) and can_reuse(snps_tranche_path):
            snp_tranche_fragments.append(get_batch().read_input(str(snps_tranche_path)))
            continue

        snps_recal_j = snps_recalibrator_scattered(
            get_batch(),
            siteonly_vcf=vcf_resources[idx],
            model_file=snp_model_in_batch,
            hapmap_resource_vcf=resources['hapmap'],
            omni_resource_vcf=resources['omni'],
            one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=50,
            use_as_annotations=True,
            job_attrs={'part': f'{idx + 1}/{fragment_count}'},
        )
        scatter_jobs.append(snps_recal_j)

        # write the results out to GCP
        get_batch().write_output(snps_recal_j.recalibration, str(snps_recal_paths[idx]))
        get_batch().write_output(snps_recal_j.recalibration_idx, str(snps_recal_paths[idx]) + '.idx')
        get_batch().write_output(snps_recal_j.tranches, str(snps_tranches_paths[idx]))

    gather_tranches_j = snps_gather_tranches_job(
        get_batch(),
        tranches=snp_tranche_fragments,
        disk_size=100,
    )
    gather_tranches_j.depends_on(*scatter_jobs)
    scatter_jobs.append(gather_tranches_j)
    get_batch().write_output(gather_tranches_j.out_tranches, output_path)
    return scatter_jobs
