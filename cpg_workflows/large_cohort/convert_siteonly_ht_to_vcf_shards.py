import logging

import hail as hl

from cpg_utils.config import config_retrieve

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def run(
    vqsr_siteonly_ht_path: str,
    repartitioned_dir: str,
    out_vcf_dir: str,
):
    vqsr_siteonly_ht = hl.read_table(str(vqsr_siteonly_ht_path))
    logger.info(f'Counting variants in VQSR site-only HT: {vqsr_siteonly_ht_path}')
    num_variants = vqsr_siteonly_ht.count()
    logger.info(f'Number of variants in VQSR site-only HT: {num_variants}')

    # will most likely need to repartition.
    # want roughly 500.000 variants per partition or if table has less than 500.000 variants, one partition
    variants_per_vep: int = config_retrieve(['workflow', 'variants_per_vep'], 500_000)
    num_partitions = max(1, num_variants // variants_per_vep)
    logger.info(f'Repartitioning VQSR site-only HT to {num_partitions} partitions')

    vqsr_siteonly_ht = vqsr_siteonly_ht.repartition(num_partitions, shuffle=True)
    vqsr_siteonly_ht = vqsr_siteonly_ht.checkpoint(
        f'{repartitioned_dir}/vqsr_siteonly_ht_repartitioned.ht',
        overwrite=True,
    )
    logger.info(f'Repartitioned VQSR site-only HT to {num_partitions} partitions')

    # Export to VCF shards
    hl.export_vcf(
        vqsr_siteonly_ht,
        out_vcf_dir,
        tabix=True,
        parallel='header_per_shard',
    )
