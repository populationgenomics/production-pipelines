"""
Create Hail Batch jobs to create and apply a VQSR models.
"""

from typing import List, Optional
import logging
import hailtop.batch as hb
from cpg_utils.hail_batch import reference_path, image_path

from cpg_pipes.hb.resources import STANDARD, HIGHMEM
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes import utils
from cpg_pipes.jobs.picard import get_intervals
from cpg_pipes.jobs.vcf import gather_vcfs
from cpg_pipes.hb.command import wrap_command

logger = logging.getLogger(__file__)


# VQSR - when applying model - targets indel_filter_level and snp_filter_level
# sensitivities. The tool matches them internally to a VQSLOD score cutoff
# based on the model's estimated sensitivity to a set of true variants.
SNP_HARD_FILTER_LEVEL = 99.7
INDEL_HARD_FILTER_LEVEL = 99.0

STANDARD_FEATURES = [
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'FS',
    'SOR',
    # 'DP',
]
SNP_STANDARD_FEATURES = STANDARD_FEATURES + ['MQ']
INDEL_STANDARD_FEATURES = STANDARD_FEATURES

ALLELE_SPECIFIC_FEATURES = [
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    # Not using depth for the following reaasons:
    # 1. The Broad pipelines don't use it;
    # 2. -G AS_StandardAnnotation flag to GenotypeGVCFs doesn't include it;
    # 3. For exomes, depth is an irrelevant feature and should be skipped.
    # 'AS_VarDP',
]
SNP_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES + ['AS_MQ']
INDEL_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES

SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]

DEFAULT_INTERVALS_NUM = 50


def make_vqsr_jobs(
    b: hb.Batch,
    input_vcf_or_mt_path: Path,
    tmp_prefix: Path,
    gvcf_count: int,
    meta_ht_path: Path | None = None,
    hard_filter_ht_path: Path | None = None,
    out_path: Path | None = None,
    use_as_annotations: bool = True,
    overwrite: bool = False,
    scatter_count: int = DEFAULT_INTERVALS_NUM,
    intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Add jobs that perform the allele-specific VQSR variant QC

    @param b: Batch object to add jobs to
    @param input_vcf_or_mt_path: path to a site-only VCF, or a matrix table
    @param meta_ht_path: if input_vcf_or_mt_path is a matrix table, this table will
           be used as a source of annotations for that matrix table, i.e.
           to filter out samples flagged as `meta.related`
    @param hard_filter_ht_path: if input_vcf_or_mt_path is a matrix table, this table
           will be used as a list of samples to hard filter out
    @param tmp_prefix: bucket for intermediate files
    @param gvcf_count: number of input samples. Can't read from combined_mt_path as it
           might not be yet genereated the point of Batch job submission
    @param scatter_count: number of intervals to parallelise SNP model creation
    @param intervals_path: path to specific interval list
    @param out_path: path to write final recalibrated VCF to
    @param use_as_annotations: use allele-specific annotation for VQSR
    @param overwrite: whether to not reuse intermediate files
    @param job_attrs: default job attributes
    @return: a final Job, and a path to the VCF with VQSR annotations
    """
    resources = {
        key: b.read_input_group(
            base=str(reference_path(f'broad/{key}_vcf')),
            index=str(reference_path(f'broad/{key}_vcf_index')),
        )
        for key in [
            'dbsnp',
            'hapmap',
            'omni',
            'one_thousand_genomes',
            'mills',
            'axiom_poly',
        ]
    }

    is_small_callset = gvcf_count < 1000
    # For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = gvcf_count >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    # To fit only a site-only VCF
    small_disk = 50 if is_small_callset else (100 if not is_huge_callset else 200)
    # To fit a joint-called VCF
    medium_disk = 100 if is_small_callset else (200 if not is_huge_callset else 500)
    huge_disk = 200 if is_small_callset else (500 if not is_huge_callset else 2000)

    jobs: list[Job] = []

    if input_vcf_or_mt_path.name.endswith('.mt'):
        # Importing dynamically to make sure $CPG_DATASET_GCP_PROJECT is set.
        from analysis_runner import dataproc

        assert meta_ht_path
        assert hard_filter_ht_path
        job_name = 'VQSR: MT to site-only VCF'
        combined_vcf_path = tmp_prefix / 'input.vcf.gz'
        if not utils.can_reuse(combined_vcf_path, overwrite):
            mt_to_vcf_job = dataproc.hail_dataproc_job(
                b,
                f'cpg_pipes/dataproc_scripts/mt_to_siteonlyvcf.py --overwrite '
                f'--mt {input_vcf_or_mt_path} '
                f'--meta-ht {meta_ht_path} '
                f'--hard-filtered-samples-ht {hard_filter_ht_path} '
                f'-o {combined_vcf_path} ',
                max_age='8h',
                packages=utils.DATAPROC_PACKAGES,
                num_secondary_workers=scatter_count,
                # hl.export_vcf() uses non-preemptible workers' disk to merge VCF files.
                # 10 samples take 2.3G, 400 samples take 60G, which roughly matches
                # `huge_disk` (also used in the AS-VQSR VCF-gather job)
                worker_boot_disk_size=huge_disk,
                job_name=job_name,
            )
        else:
            mt_to_vcf_job = b.new_job(f'{job_name} [reuse]', job_attrs)
        jobs.append(mt_to_vcf_job)
        tabix_job = add_tabix_job(
            b,
            vcf_path=combined_vcf_path,
            disk_size=medium_disk,
            job_attrs=job_attrs,
        )
        tabix_job.depends_on(mt_to_vcf_job)
        siteonly_vcf = tabix_job.combined_vcf
    else:
        siteonly_vcf = b.read_input_group(
            **{
                'vcf.gz': str(input_vcf_or_mt_path),
                'vcf.gz.tbi': str(input_vcf_or_mt_path) + '.tbi',
            }
        )

    indels_variant_recalibrator_job = add_indels_variant_recalibrator_job(
        b=b,
        sites_only_variant_filtered_vcf=siteonly_vcf,
        mills_resource_vcf=resources['mills'],
        axiom_poly_resource_vcf=resources['axiom_poly'],
        dbsnp_resource_vcf=resources['dbsnp'],
        disk_size=small_disk,
        use_as_annotations=use_as_annotations,
        is_small_callset=is_small_callset,
        job_attrs=job_attrs,
    )
    jobs.append(indels_variant_recalibrator_job)

    indels_recalibration = indels_variant_recalibrator_job.recalibration
    indels_tranches = indels_variant_recalibrator_job.tranches

    snp_max_gaussians = 6
    if is_small_callset:
        snp_max_gaussians = 4
    elif is_huge_callset:
        snp_max_gaussians = 8

    if scatter_count > 1:
        intervals_j, intervals = get_intervals(
            b=b,
            scatter_count=scatter_count,
            intervals_path=intervals_path,
            job_attrs=job_attrs,
            output_prefix=tmp_prefix,
        )
        jobs.append(intervals_j)

        # Run SNP recalibrator in a scattered mode
        model_j = add_snps_variant_recalibrator_create_model_step(
            b,
            sites_only_variant_filtered_vcf=siteonly_vcf,
            hapmap_resource_vcf=resources['hapmap'],
            omni_resource_vcf=resources['omni'],
            one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            max_gaussians=snp_max_gaussians,
            job_attrs=job_attrs,
        )
        jobs.append(model_j)

        snps_recalibrator_jobs = [
            add_snps_variant_recalibrator_scattered_step(
                b,
                sites_only_vcf=siteonly_vcf,
                interval=intervals[idx],
                model_file=model_j.model_file,
                hapmap_resource_vcf=resources['hapmap'],
                omni_resource_vcf=resources['omni'],
                one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
                dbsnp_resource_vcf=resources['dbsnp'],
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                max_gaussians=snp_max_gaussians,
                is_small_callset=is_small_callset,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            for idx in range(scatter_count)
        ]
        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = add_snps_gather_tranches_step(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
        ).out_tranches

        scattered_vcfs = [
            add_apply_recalibration_step(
                b,
                input_vcf=siteonly_vcf,
                interval=intervals[idx],
                indels_recalibration=indels_recalibration,
                indels_tranches=indels_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                disk_size=huge_disk,
                use_as_annotations=use_as_annotations,
                snp_filter_level=SNP_HARD_FILTER_LEVEL,
                indel_filter_level=INDEL_HARD_FILTER_LEVEL,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            ).output_vcf
            for idx in range(scatter_count)
        ]
        recalibrated_gathered_vcf_jobs, recalibrated_gathered_vcf = gather_vcfs(
            b=b,
            input_vcfs=[v['vcf.gz'] for v in scattered_vcfs],
            overwrite=overwrite,
            out_vcf_path=out_path,
            site_only=True,
            gvcf_count=gvcf_count,
            job_attrs=job_attrs,
        )
        for j in recalibrated_gathered_vcf_jobs:
            j.name = f'VQSR: {j.name}'
        jobs.extend(recalibrated_gathered_vcf_jobs)

    else:
        snps_recalibrator_job = add_snps_variant_recalibrator_step(
            b,
            sites_only_variant_filtered_vcf=siteonly_vcf,
            hapmap_resource_vcf=resources['hapmap'],
            omni_resource_vcf=resources['omni'],
            one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            max_gaussians=snp_max_gaussians,
            is_small_callset=is_small_callset,
            job_attrs=job_attrs,
        )
        jobs.append(snps_recalibrator_job)

        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        recalibrated_gathered_vcf_j = add_apply_recalibration_step(
            b,
            input_vcf=siteonly_vcf,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            indel_filter_level=SNP_HARD_FILTER_LEVEL,
            snp_filter_level=INDEL_HARD_FILTER_LEVEL,
            output_vcf_path=out_path,
            job_attrs=job_attrs,
        )
        jobs.append(recalibrated_gathered_vcf_j)

    return jobs


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'VQSR: MakeSitesOnlyVcf'
    if utils.can_reuse(output_vcf_path, overwrite):
        return b.new_job(job_name + ' [reuse]')

    j = b.new_job(job_name)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, mem_gb=8, storage_gb=disk_size)
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        wrap_command(
            f"""\
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    MakeSitesOnlyVcf \\
    -I {input_vcf['vcf.gz']} \\
    -O {j.output_vcf['vcf.gz']} \\
    --CREATE_INDEX
    """
        )
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def add_tabix_job(
    b: hb.Batch,
    vcf_path: Path,
    disk_size: int,
    job_attrs: dict | None = None,
) -> Job:
    """
    Regzip and tabix the combined VCF (for some reason the one produced by
    `mt_to_vcf.py` is not block-gzipped).
    """
    j = b.new_job('VQSR: Tabix', job_attrs)
    j.image(image_path('bcftools'))
    STANDARD.set_resources(j, mem_gb=8, storage_gb=disk_size)
    j.declare_resource_group(
        combined_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    vcf_inp = b.read_input(str(vcf_path))
    j.command(
        wrap_command(
            f"""\
    gunzip {vcf_inp} -c | bgzip -c > {j.combined_vcf['vcf.gz']}
    tabix -p vcf {j.combined_vcf['vcf.gz']}
    """
        )
    )
    return j


def add_sites_only_gather_vcf_step(
    b: hb.Batch,
    input_vcfs: list[hb.ResourceGroup],
    disk_size: int,
    job_attrs: dict | None = None,
) -> Job:
    """
    Gathers VCF files from scattered operations into a single VCF file

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    j = b.new_job('VQSR: SitesOnlyGatherVcf', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, mem_gb=8, storage_gb=disk_size)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(
        wrap_command(
            f"""\
    # --ignore-safety-checks makes a big performance difference so we include it in
    # our invocation. This argument disables expensive checks that the file headers
    # contain the same set of genotyped samples and that files are in order by position
    # of first record.
    gatk --java-options -Xms{res.get_java_mem_mb()}m GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    """
        )
    )
    return j


def add_indels_variant_recalibrator_job(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    mills_resource_vcf: hb.ResourceGroup,
    axiom_poly_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    max_gaussians: int = 4,
    is_small_callset: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run VariantRecalibrator to calculate VQSLOD tranches for indels

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value. 4 is a
    reasonable default for indels, as their number is smaller than SNPs.

    Returns: a Job object with 2 outputs: j.recalibration (ResourceGroup), j.tranches.
    """
    j = b.new_job('VQSR: IndelsVariantRecalibrator', job_attrs)
    j.image(image_path('gatk'))

    # We run it for the entire dataset in one job, so can take an entire instance.
    instance_fraction = 1
    # however, for smaller datasets we take a standard instance, and for larger
    # ones we take a highmem instance
    if is_small_callset:
        res = STANDARD.set_resources(
            j, fraction=instance_fraction, storage_gb=disk_size
        )
    else:
        res = HIGHMEM.set_resources(j, fraction=instance_fraction, storage_gb=disk_size)

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join(
        [f'-tranche {v}' for v in INDEL_RECALIBRATION_TRANCHE_VALUES]
    )
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                INDEL_ALLELE_SPECIFIC_FEATURES
                if use_as_annotations
                else INDEL_STANDARD_FEATURES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode INDEL \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --max-gaussians {max_gaussians} \\
      -resource:mills,known=false,training=true,truth=true,prior=12 {mills_resource_vcf.base} \\
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {axiom_poly_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 {dbsnp_resource_vcf.base}
      """
    )
    return j


def add_snps_variant_recalibrator_create_model_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    is_small_callset: bool = False,
    is_huge_callset: bool = False,
    max_gaussians: int = 4,
    job_attrs: dict | None = None,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 1 output j.model
    The latter is useful to produce the optional tranche plot.
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorCreateModel', job_attrs)
    j.image(image_path('gatk'))

    # We run it for the entire dataset in one job, so can take an entire instance.
    instance_fraction = 1
    # however, for smaller datasets we take a standard instance, and for larger
    # ones we take a highmem instance
    if is_small_callset:
        res = STANDARD.set_resources(
            j, fraction=instance_fraction, storage_gb=disk_size
        )
    else:
        res = HIGHMEM.set_resources(j, fraction=instance_fraction, storage_gb=disk_size)

    downsample_factor = 75 if is_huge_callset else 10

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_ALLELE_SPECIFIC_FEATURES
                if use_as_annotations
                else SNP_STANDARD_FEATURES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --sample-every-Nth-variant {downsample_factor} \\
      --output-model {j.model_file} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}
      """
    )
    return j


def add_snps_variant_recalibrator_scattered_step(
    b: hb.Batch,
    sites_only_vcf: hb.ResourceGroup,
    model_file: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    interval: hb.Resource | None = None,
    max_gaussians: int = 4,
    is_small_callset: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    Returns: a Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches
    """
    j = b.new_job('VQSR: SNPsVariantRecalibratorScattered', job_attrs)
    j.image(image_path('gatk'))

    if is_small_callset:
        res = STANDARD.set_resources(j, ncpu=4, storage_gb=disk_size)
    else:
        res = STANDARD.set_resources(j, ncpu=8, storage_gb=disk_size)

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_ALLELE_SPECIFIC_FEATURES
                if use_as_annotations
                else SNP_STANDARD_FEATURES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    MODEL_REPORT={model_file}

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      VariantRecalibrator \\
      -V {sites_only_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {f'-L {interval} ' if interval else ''} \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --input-model {model_file} --output-tranches-for-scatter \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}"""
    )
    return j


def add_snps_variant_recalibrator_step(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    max_gaussians: int = 4,
    is_small_callset: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    """
    j = b.new_job('VQSR: SNPsVariantRecalibrator', job_attrs)
    j.image(image_path('gatk'))

    # We run it for the entire dataset in one job, so can take an entire instance.
    instance_fraction = 1
    # however, for smaller datasets we take a standard instance, and for larger
    # ones we take a highmem instance
    if is_small_callset:
        res = STANDARD.set_resources(
            j, fraction=instance_fraction, storage_gb=disk_size
        )
    else:
        res = HIGHMEM.set_resources(j, fraction=instance_fraction, storage_gb=disk_size)

    j.storage(f'{disk_size}G')

    j.declare_resource_group(recalibration={'index': '{root}.idx', 'base': '{root}'})

    tranche_cmdl = ' '.join([f'-tranche {v}' for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = ' '.join(
        [
            f'-an {v}'
            for v in (
                SNP_ALLELE_SPECIFIC_FEATURES
                if use_as_annotations
                else SNP_STANDARD_FEATURES
            )
        ]
    )
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      VariantRecalibrator \\
      -V {sites_only_variant_filtered_vcf['vcf.gz']} \\
      -O {j.recalibration} \\
      --tranches-file {j.tranches} \\
      --trust-all-polymorphic \\
      {tranche_cmdl} \\
      {an_cmdl} \\
      -mode SNP \\
      {"--use-allele-specific-annotations " if use_as_annotations else ""} \\
      --max-gaussians {max_gaussians} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15 {hapmap_resource_vcf.base} \\
      -resource:omni,known=false,training=true,truth=true,prior=12 {omni_resource_vcf.base} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10 {one_thousand_genomes_resource_vcf.base} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base} 
      """
    )
    return j


def add_snps_gather_tranches_step(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    disk_size: int,
    job_attrs: dict | None = None,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibrartion process
    is borken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    Returns: a Job object with one output j.out_tranches
    """
    j = b.new_job('VQSR: SNPGatherTranches', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)

    inputs_cmdl = ' '.join([f'--input {t}' for t in tranches])
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      GatherTranches \\
      --mode SNP \\
      {inputs_cmdl} \\
      --output {j.out_tranches}"""
    )
    return j


def add_apply_recalibration_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    indels_recalibration: hb.ResourceGroup,
    indels_tranches: hb.ResourceFile,
    snps_recalibration: hb.ResourceGroup,
    snps_tranches: hb.ResourceFile,
    disk_size: int,
    use_as_annotations: bool,
    indel_filter_level: float,
    snp_filter_level: float,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.
    Runs ApplyVQSR twice to apply first indel, then SNP recalibrations.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truthset. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    Returns: a Job object with one ResourceGroup output j.output_vcf, corresponding
    to a VCF with tranche annotated in the FILTER field
    """
    j = b.new_job('VQSR: ApplyVQSR', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(
        f"""set -euo pipefail

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    TMP_DIR=$(dirname {j.output_vcf['vcf.gz']})/tmp
    mkdir $TMP_DIR

    TMP_INDEL_RECALIBRATED=/io/batch/tmp.indel.recalibrated.vcf.gz
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      ApplyVQSR \\
      --tmp-dir $TMP_DIR \\
      -O $TMP_INDEL_RECALIBRATED \\
      -V {input_vcf['vcf.gz']} \\
      --recal-file {indels_recalibration} \\
      --tranches-file {indels_tranches} \\
      --truth-sensitivity-filter-level {indel_filter_level} \\
      --create-output-variant-index true \\
      {f'-L {interval} ' if interval else ''} \\
      {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
      -mode INDEL
      
    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    rm {input_vcf['vcf.gz']} {indels_recalibration} {indels_tranches}
    rm -rf $TMP_DIR
    mkdir $TMP_DIR

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      ApplyVQSR \\
      -O {j.output_vcf['vcf.gz']} \\
      -V $TMP_INDEL_RECALIBRATED \\
      --recal-file {snps_recalibration} \\
      --tranches-file {snps_tranches} \\
      --truth-sensitivity-filter-level {snp_filter_level} \\
      --create-output-variant-index true \\
      {f'-L {interval} ' if interval else ''} \\
      {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
      -mode SNP

    df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
      """
    )

    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j


def add_collect_metrics_sharded_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    interval_list: hb.ResourceFile,
    ref_dict: hb.ResourceFile,
    disk_size: int,
    job_attrs: dict | None = None,
):
    """
    Run CollectVariantCallingMetrics for site-level evaluation.

    This method produces detailed and summary metrics report files. The summary metrics
    provide cohort-level variant metrics and the detailed metrics segment variant
    metrics for each sample in the callset. The detail metrics give the same metrics
    as the summary metrics for the samples plus several additional metrics.

    These are explained in detail at
    https://broadinstitute.github.io/picard/picard-metric-definitions.html.

    Returns: a `Job` object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles
    """
    j = b.new_job('VQSR: CollectMetricsSharded', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      CollectVariantCallingMetrics \\
      --INPUT {input_vcf['vcf.gz']} \\
      --DBSNP {dbsnp_vcf.base} \\
      --SEQUENCE_DICTIONARY {ref_dict} \\
      --OUTPUT {j.metrics} \\
      --THREAD_COUNT 8 \\
      --TARGET_INTERVALS {interval_list}"""
    )
    return j


def _add_variant_eval_step(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    ref_fasta: hb.ResourceGroup,
    dbsnp_vcf: hb.ResourceGroup,
    disk_size: int,
    output_path: str = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run VariantEval for site-level evaluation.
    Saves the QC to `output_path` bucket
    """
    j = b.new_job('VQSR: VariantEval', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)

    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      VariantEval \\
      --eval {input_vcf['vcf.gz']} \\
      -R {ref_fasta.base} \\
      -D {dbsnp_vcf.base} \\
      --output {j.output}"""
    )
    if output_path:
        b.write_output(j.output, output_path)
    return j


def add_gather_variant_calling_metrics_step(
    b: hb.Batch,
    input_details: List[hb.ResourceGroup],
    input_summaries: List[hb.ResourceGroup],
    disk_size: int,
    output_path_prefix: str = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Combines metrics from multiple CollectVariantCallingMetrics runs.

    Returns: a `Job` object with a single ResourceGroup output j.metrics, with
    j.metrics.detail_metrics and j.metrics.summary_metrics ResourceFiles

    Saves the QC results to a bucket with the `output_path_prefix` prefix
    """
    j = b.new_job('VQSR: GatherVariantCallingMetrics', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)
    j.declare_resource_group(
        metrics={
            'detail_metrics': '{root}.variant_calling_detail_metrics',
            'summary_metrics': '{root}.variant_calling_summary_metrics',
        }
    )

    input_cmdl = ' '.join(f'--INPUT {f} ' for f in input_details + input_summaries)
    j.command(
        f"""set -euo pipefail

    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
      AccumulateVariantCallingMetrics \\
      {input_cmdl} \\
      --OUTPUT {j.metrics}"""
    )
    if output_path_prefix:
        b.write_output(j.metrics, output_path_prefix)
    return j
