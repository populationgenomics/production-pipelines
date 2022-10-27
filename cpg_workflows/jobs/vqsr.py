"""
Create Hail Batch jobs to create and apply AS-VQSR models.

Parameters are borrowed from WARP:
WGS VQSR: https://github.com/broadinstitute/warp/blob/79261cde9bd06bb6b1d4a83d75dc54f734541fec/pipelines/broad/dna_seq/germline/joint_genotyping/wgs/JointGenotyping.inputs.json#L29-L35 (there is no direct example config for WGS AS-VQSR, but adjusted correspondingly)
Exome AS-VQSR: https://github.com/broadinstitute/warp/blob/79261cde9bd06bb6b1d4a83d75dc54f734541fec/pipelines/broad/dna_seq/germline/joint_genotyping/exome/JointGenotyping.inputs.json#L8-L11
Note that there is no example settings config for WGS AS-VQSR, so we construct it 
from WGS VQSR and Exome AS-VQSR settings.
"""

from typing import List

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, image_path, command
from cpg_workflows.resources import STANDARD, HIGHMEM
from cpg_workflows.jobs.picard import get_intervals
from cpg_workflows.jobs.vcf import gather_vcfs, subset_vcf

STANDARD_FEATURES = [
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'FS',
    'SOR',
]
SNP_STANDARD_FEATURES = STANDARD_FEATURES + ['MQ']
INDEL_STANDARD_FEATURES = STANDARD_FEATURES

ALLELE_SPECIFIC_FEATURES = [
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    # Not using depth for the following reasons:
    # 1. The Broad pipelines don't use it;
    # 2. -G AS_StandardAnnotation flag to GenotypeGVCFs doesn't include it;
    # 3. For exomes, depth is an irrelevant feature and should be skipped:
    # 'AS_VarDP'
    # Note that for consistency, we also skip it for WGS.
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


def make_vqsr_jobs(
    b: hb.Batch,
    input_siteonly_vcf_path: Path,
    tmp_prefix: Path,
    gvcf_count: int,
    out_path: Path | None = None,
    use_as_annotations: bool = True,
    overwrite: bool = False,
    intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Add jobs that perform the allele-specific VQSR variant QC

    @param b: Batch object to add jobs to
    @param input_siteonly_vcf_path: path to a site-only VCF
    @param tmp_prefix: bucket for intermediate files
    @param gvcf_count: number of input samples. Can't read from combined_mt_path as it
           might not be yet generated the point of Batch job submission
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

    # When applying model, VQSR targets indel_filter_level and snp_filter_level
    # sensitivities. The tool matches them internally to a VQSLOD score cutoff
    # based on the model's estimated sensitivity to a set of true variants.
    snp_filter_level = get_config()['vqsr']['snp_filter_level']
    indel_filter_level = get_config()['vqsr']['indel_filter_level']

    is_small_callset = gvcf_count < 1000
    # For small callsets, we gather the VCF shards and collect QC metrics directly.
    # For anything larger, we need to keep the VCF sharded and gather metrics
    # collected from them.
    is_huge_callset = gvcf_count >= 100000
    # For huge callsets, we allocate more memory for the SNPs Create Model step

    # To fit only a site-only VCF
    small_disk = 50 if is_small_callset else (100 if not is_huge_callset else 200)
    # To fit a joint-called VCF
    huge_disk = 200 if is_small_callset else (500 if not is_huge_callset else 2000)

    scatter_count = 50
    if gvcf_count > 300:
        scatter_count = 100
    if gvcf_count > 1000:
        scatter_count = 200

    jobs: list[Job] = []

    siteonly_vcf = b.read_input_group(
        **{
            'vcf.gz': str(input_siteonly_vcf_path),
            'vcf.gz.tbi': str(input_siteonly_vcf_path) + '.tbi',
        }
    )

    indel_vcf_j = subset_vcf(
        b=b,
        vcf=siteonly_vcf,
        variant_types=['INDEL', 'MNP', 'MIXED'],
        job_attrs=job_attrs,
    )
    jobs.append(indel_vcf_j)
    indel_vcf_j.name = f'VQSR: {indel_vcf_j.name}'
    indel_vcf = indel_vcf_j.output_vcf
    assert isinstance(indel_vcf, hb.ResourceGroup)
    indel_recalibrator_j = indel_recalibrator_job(
        b=b,
        siteonly_vcf=indel_vcf,
        mills_resource_vcf=resources['mills'],
        axiom_poly_resource_vcf=resources['axiom_poly'],
        dbsnp_resource_vcf=resources['dbsnp'],
        disk_size=small_disk,
        use_as_annotations=use_as_annotations,
        is_small_callset=is_small_callset,
        job_attrs=job_attrs,
    )
    # To make type checkers happy:
    assert isinstance(indel_recalibrator_j.recalibration, hb.ResourceGroup)
    assert isinstance(indel_recalibrator_j.tranches, hb.ResourceFile)
    jobs.append(indel_recalibrator_j)

    if scatter_count > 1:
        # Run SNP recalibrator in a scattered mode
        model_j = snps_recalibrator_create_model_job(
            b,
            siteonly_vcf=siteonly_vcf,
            hapmap_resource_vcf=resources['hapmap'],
            omni_resource_vcf=resources['omni'],
            one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            is_small_callset=is_small_callset,
            is_huge_callset=is_huge_callset,
            job_attrs=job_attrs,
        )
        jobs.append(model_j)
        assert isinstance(model_j.model_file, hb.ResourceFile)

        intervals_j, intervals = get_intervals(
            b=b,
            scatter_count=scatter_count,
            source_intervals_path=intervals_path,
            job_attrs=job_attrs,
            output_prefix=tmp_prefix / f'intervals_{scatter_count}',
        )
        if intervals_j:
            jobs.append(intervals_j)

        snps_interval_vcfs = []
        for idx in range(scatter_count):
            snp_subset_j = subset_vcf(
                b=b,
                vcf=siteonly_vcf,
                interval=intervals[idx],
                variant_types=['SNP'],
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            snp_subset_j.name = f'VQSR: {snp_subset_j.name}'
            jobs.append(snp_subset_j)
            snp_subset_vcf = snp_subset_j.output_vcf
            snps_interval_vcfs.append(snp_subset_vcf)

        snps_recal_jobs = []
        for idx in range(scatter_count):
            snps_interval_vcf = snps_interval_vcfs[idx]
            assert isinstance(snps_interval_vcf, hb.ResourceGroup)
            snps_recal_j = snps_recalibrator_scattered(
                b,
                siteonly_vcf=snps_interval_vcf,
                interval=intervals[idx],
                model_file=model_j.model_file,
                hapmap_resource_vcf=resources['hapmap'],
                omni_resource_vcf=resources['omni'],
                one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
                dbsnp_resource_vcf=resources['dbsnp'],
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                is_small_callset=is_small_callset,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            snps_recal_jobs.append(snps_recal_j)

        snps_recalibrations = [j.recalibration for j in snps_recal_jobs]
        snps_tranches = []
        for snps_recal_j in snps_recal_jobs:
            assert isinstance(snps_recal_j.tranches, hb.ResourceFile)
            snps_tranches.append(snps_recal_j.tranches)
        snp_gathered_tranches = snps_gather_tranches_job(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
            job_attrs=job_attrs,
        ).out_tranches

        assert isinstance(snp_gathered_tranches, hb.ResourceFile)
        snps_interval_snp_applied_vcfs = []
        for idx in range(scatter_count):
            snp_recalibration = snps_recalibrations[idx]
            assert isinstance(snp_recalibration, hb.ResourceGroup)
            snps_interval_vcf = snps_interval_vcfs[idx]
            assert isinstance(snps_interval_vcf, hb.ResourceGroup)
            applied_snps_vcf = apply_recalibration_snps(
                b,
                input_vcf=snps_interval_vcf,
                interval=intervals[idx],
                recalibration=snp_recalibration,
                tranches=snp_gathered_tranches,
                disk_size=huge_disk,
                use_as_annotations=use_as_annotations,
                filter_level=snp_filter_level,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            ).output_vcf
            assert isinstance(applied_snps_vcf, hb.ResourceGroup)
            snps_interval_snp_applied_vcfs.append(applied_snps_vcf['vcf.gz'])

        snps_applied_gathered_j, snps_applied_gathered_vcf = gather_vcfs(
            b=b,
            input_vcfs=snps_interval_snp_applied_vcfs + [indel_vcf['vcf.gz']],
            overwrite=overwrite,
            site_only=True,
            gvcf_count=gvcf_count,
            job_attrs=job_attrs,
        )
        if snps_applied_gathered_j:
            snps_applied_gathered_j.name = f'VQSR: {snps_applied_gathered_j.name}'
            jobs.append(snps_applied_gathered_j)

        apply_indel_j = apply_recalibration_indels(
            b,
            input_vcf=snps_applied_gathered_vcf,
            recalibration=indel_recalibrator_j.recalibration,
            tranches=indel_recalibrator_j.tranches,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            filter_level=indel_filter_level,
            job_attrs=job_attrs,
        )
        jobs.append(apply_indel_j)
        final_vcf = apply_indel_j.output_vcf

    else:
        snps_recal_j = snps_recalibrator_job(
            b,
            sites_only_variant_filtered_vcf=siteonly_vcf,
            hapmap_resource_vcf=resources['hapmap'],
            omni_resource_vcf=resources['omni'],
            one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            is_small_callset=is_small_callset,
            job_attrs=job_attrs,
        )
        jobs.append(snps_recal_j)

        assert isinstance(snps_recal_j.recalibration, hb.ResourceGroup)
        assert isinstance(snps_recal_j.tranches, hb.ResourceFile)
        assert isinstance(snps_recal_j.recalibration, hb.ResourceGroup)
        assert isinstance(snps_recal_j.tranches, hb.ResourceFile)

        apply_recal_snps_j = apply_recalibration_snps(
            b,
            input_vcf=siteonly_vcf,
            recalibration=snps_recal_j.recalibration,
            tranches=snps_recal_j.tranches,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            filter_level=snp_filter_level,
            job_attrs=job_attrs,
        )
        assert isinstance(apply_recal_snps_j.output_vcf, hb.ResourceGroup)
        apply_indel_j = apply_recalibration_indels(
            b,
            input_vcf=apply_recal_snps_j.output_vcf,
            recalibration=indel_recalibrator_j.recalibration,
            tranches=indel_recalibrator_j.tranches,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            filter_level=indel_filter_level,
            job_attrs=job_attrs,
        )
        jobs.append(apply_indel_j)
        final_vcf = apply_indel_j.output_vcf

    if out_path:
        b.write_output(final_vcf, str(out_path).replace('.vcf.gz', ''))
    return jobs


def add_tabix_job(
    b: hb.Batch,
    vcf_path: Path,
    disk_size: int,
    job_attrs: dict | None = None,
) -> Job:
    """
    Re-gzip and tabix the combined VCF (for some reason the one produced by
    `mt_to_vcf.py` is not block-gzipped).
    """
    job_attrs = (job_attrs or {}) | {'tool': 'tabix'}
    j = b.new_job('VQSR: Tabix', job_attrs)
    j.image(image_path('bcftools'))
    STANDARD.set_resources(j, mem_gb=8, storage_gb=disk_size)
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    vcf_inp = b.read_input(str(vcf_path))
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    j.command(
        command(
            f"""\
    gunzip {vcf_inp} -c | bgzip -c > {j.output_vcf['vcf.gz']}
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """
        )
    )
    return j


def indel_recalibrator_job(
    b: hb.Batch,
    siteonly_vcf: hb.ResourceGroup,
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
    job_attrs = (job_attrs or {}) | {'tool': 'gatk VariantRecalibrator'}
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
      -V {siteonly_vcf['vcf.gz']} \\
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


def snps_recalibrator_create_model_job(
    b: hb.Batch,
    siteonly_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    is_small_callset: bool = False,
    is_huge_callset: bool = False,
    max_gaussians: int = 6,
    job_attrs: dict | None = None,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants
    and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:
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
    job_attrs = (job_attrs or {}) | {'tool': 'gatk VariantRecalibrator'}
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
      -V {siteonly_vcf['vcf.gz']} \\
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


def snps_recalibrator_scattered(
    b: hb.Batch,
    siteonly_vcf: hb.ResourceGroup,
    model_file: hb.ResourceFile,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    interval: hb.Resource | None = None,
    max_gaussians: int = 6,
    is_small_callset: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply
    the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:
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
    job_attrs = (job_attrs or {}) | {'tool': 'gatk VariantRecalibrator'}
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
      -V {siteonly_vcf['vcf.gz']} \\
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


def snps_recalibrator_job(
    b: hb.Batch,
    sites_only_variant_filtered_vcf: hb.ResourceGroup,
    hapmap_resource_vcf: hb.ResourceGroup,
    omni_resource_vcf: hb.ResourceGroup,
    one_thousand_genomes_resource_vcf: hb.ResourceGroup,
    dbsnp_resource_vcf: hb.ResourceGroup,
    disk_size: int,
    use_as_annotations: bool,
    max_gaussians: int = 6,
    is_small_callset: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Recalibrate SNPs in one run (alternative to scatter-gather approach)
    """
    job_attrs = (job_attrs or {}) | {'tool': 'gatk VariantRecalibrator'}
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


def snps_gather_tranches_job(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    disk_size: int,
    job_attrs: dict | None = None,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval
    tranches outputs.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:
    1. Run the recalibrator with the following additional arguments:
       --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
    2. Apply the resulting model to each genomic interval with, running the recalibrator
       with the same base parameters, plus:
       --input-model <model-file> --output-tranches-for-scatter
    3. Collate the resulting per-interval tranches with GatherTranches

    Returns: a Job object with one output j.out_tranches
    """
    job_attrs = (job_attrs or {}) | {'tool': 'gatk GatherTranches'}
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


def apply_recalibration_snps(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    recalibration: hb.ResourceGroup,
    tranches: hb.ResourceFile,
    disk_size: int,
    use_as_annotations: bool,
    filter_level: float,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truth-set. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    Returns: a Job object with one ResourceGroup output j.output_vcf, corresponding
    to a VCF with tranche annotated in the FILTER field
    """
    job_attrs = (job_attrs or {}) | {'tool': 'gatk ApplyVQSR'}
    j = b.new_job('VQSR: ApplyVQSR SNPs', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)

    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz'})
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    ApplyVQSR \\
    -O {j.output_vcf['vcf.gz']} \\
    -V {input_vcf['vcf.gz']} \\
    --recal-file {recalibration} \\
    --tranches-file {tranches} \\
    --truth-sensitivity-filter-level {filter_level} \\
    --create-output-variant-index true \\
    {f'-L {interval} ' if interval else ''} \\
    {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
    -mode SNP
    """
    j.command(command(cmd, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j


def apply_recalibration_indels(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    recalibration: hb.ResourceGroup,
    tranches: hb.ResourceFile,
    disk_size: int,
    use_as_annotations: bool,
    filter_level: float,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truth-set. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    Returns: a Job object with one ResourceGroup output j.output_vcf, corresponding
    to a VCF with tranche annotated in the FILTER field
    """
    job_attrs = (job_attrs or {}) | {'tool': 'gatk ApplyVQSR'}
    j = b.new_job('VQSR: ApplyVQSR INDEL', job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2, storage_gb=disk_size)

    j.declare_resource_group(
        output_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)

    cmd = f"""\
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    ApplyVQSR \\
    --tmp-dir $BATCH_TMPDIR \\
    -O {j.output_vcf['vcf.gz']} \\
    -V {input_vcf['vcf.gz']} \\
    --recal-file {recalibration} \\
    --tranches-file {tranches} \\
    --truth-sensitivity-filter-level {filter_level} \\
    --create-output-variant-index true \\
    {f'-L {interval} ' if interval else ''} \\
    {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
    -mode INDEL
    """
    j.command(command(cmd, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j
