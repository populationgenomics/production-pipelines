import logging

import hail as hl


def run(
    vqsr_siteonly_ht_path: str,
    tmp_dir: str,
    out_vcf_dir: str,
):
    vqsr_siteonly_ht = hl.read_table(str(vqsr_siteonly_ht_path))
    num_variants = vqsr_siteonly_ht.count()

    # will most likely need to repartition.
    # want roughly 500.000 variants per partition or if table has less than 500.000 variants, one partition
    num_partitions = max(1, num_variants // 500_000)

    vqsr_siteonly_ht = vqsr_siteonly_ht.repartition(num_partitions, shuffle=True)
    vqsr_siteonly_ht = vqsr_siteonly_ht.checkpoint(f'{tmp_dir}/vqsr_siteonly_ht_repartitioned.ht', overwrite=True)
    logging.info(f'Repartitioned VQSR site-only HT to {num_partitions} partitions')

    # Export to VCF shards
    hl.export_vcf(
        vqsr_siteonly_ht,
        out_vcf_dir,
        tabix=True,
        parallel='header_per_shard',
    )
