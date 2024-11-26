"""
Create and run jobs relating to VQSR for the RD combiner
"""

from cpg_utils.config import reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.vqsr import indel_recalibrator_job, snps_recalibrator_create_model_job


def train_vqsr_indels(sites_only_vcf: str, indel_recal: str, indel_tranches: str):
    """
    Train VQSR indels on the sites-only VCF
    Args:
        sites_only_vcf ():
        indel_recal ():
        indel_tranches ():
    """
    resources = {
        key: get_batch().read_input_group(
            base=reference_path(f'broad/{key}_vcf'),
            index=reference_path(f'broad/{key}_vcf_index'),
        )
        for key in [
            'dbsnp',
            'mills',
            'axiom_poly',
        ]
    }
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
        disk_size=100,
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
    resources = {
        key: get_batch().read_input_group(
            base=reference_path(f'broad/{key}_vcf'),
            index=reference_path(f'broad/{key}_vcf_index'),
        )
        for key in [
            'dbsnp',
            'hapmap',
            'omni',
            'one_thousand_genomes',
        ]
    }
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
        disk_size=100,
        use_as_annotations=True,
        is_small_callset=False,
    )

    get_batch().write_output(snp_recalibrator_j.model_file, snp_model)
    return snp_recalibrator_j
