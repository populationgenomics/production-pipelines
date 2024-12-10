"""
Create and run jobs relating to VQSR for the RD combiner
"""

from functools import lru_cache

from hailtop.batch.job import Job
from hailtop.batch.resource import ResourceGroup

from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.vqsr import (
    INDEL_ALLELE_SPECIFIC_FEATURES,
    INDEL_RECALIBRATION_TRANCHE_VALUES,
    SNP_ALLELE_SPECIFIC_FEATURES,
    SNP_RECALIBRATION_TRANCHE_VALUES,
)
from cpg_workflows.resources import HIGHMEM
from cpg_workflows.utils import VCF_GZ, VCF_GZ_TBI

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
