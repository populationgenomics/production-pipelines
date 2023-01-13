"""
Create Hail Batch jobs to create and apply AS-VQSR models.

Parameters are borrowed from WARP:
WGS VQSR: https://github.com/broadinstitute/warp/blob/79261cde9bd06bb6b1d4a83d75dc54f734541fec/pipelines/broad/dna_seq/germline/joint_genotyping/wgs/JointGenotyping.inputs.json#L29-L35 (there is no direct example config for WGS AS-VQSR, but adjusted correspondingly)
Exome AS-VQSR: https://github.com/broadinstitute/warp/blob/79261cde9bd06bb6b1d4a83d75dc54f734541fec/pipelines/broad/dna_seq/germline/joint_genotyping/exome/JointGenotyping.inputs.json#L8-L11
Note that there is no example settings config for WGS AS-VQSR, so we construct it 
from WGS VQSR and Exome AS-VQSR settings.
"""

from typing import List, Sequence

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, image_path, command
from cpg_workflows.resources import STANDARD, HIGHMEM, joint_calling_scatter_count
from cpg_workflows.jobs.picard import get_intervals
from cpg_workflows.jobs.vcf import gather_vcfs
from cpg_workflows.utils import can_reuse

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
    scatter_count: int | None = None,
    out_path: Path | None = None,
    use_as_annotations: bool = True,
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
    @param scatter_count: number of partitions
    @param intervals_path: path to specific interval list
    @param out_path: path to write final recalibrated VCF to
    @param use_as_annotations: use allele-specific annotation for VQSR
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

    jobs: list[Job] = []

    siteonly_vcf = b.read_input_group(
        **{
            'vcf.gz': str(input_siteonly_vcf_path),
            'vcf.gz.tbi': str(input_siteonly_vcf_path) + '.tbi',
        }
    )

    indel_recal_path = tmp_prefix / 'indel_recalibrations'
    indel_tranches_path = tmp_prefix / 'indel_tranches'
    if not can_reuse(
        [indel_recal_path, str(indel_recal_path) + '.idx', indel_tranches_path]
    ):
        indel_recalibrator_j = indel_recalibrator_job(
            b=b,
            siteonly_vcf=siteonly_vcf,
            mills_resource_vcf=resources['mills'],
            axiom_poly_resource_vcf=resources['axiom_poly'],
            dbsnp_resource_vcf=resources['dbsnp'],
            disk_size=small_disk,
            use_as_annotations=use_as_annotations,
            is_small_callset=is_small_callset,
            job_attrs=job_attrs,
        )
        jobs.append(indel_recalibrator_j)
        b.write_output(indel_recalibrator_j.recalibration, str(indel_recal_path))
        b.write_output(
            indel_recalibrator_j.recalibration_idx, str(indel_recal_path) + '.idx'
        )
        b.write_output(indel_recalibrator_j.tranches, str(indel_tranches_path))
    indel_recalibration = b.read_input(str(indel_recal_path))
    indel_recalibration_idx = b.read_input(str(indel_recal_path) + '.idx')
    indel_tranches = b.read_input(str(indel_tranches_path))

    scatter_count = scatter_count or joint_calling_scatter_count(gvcf_count)
    assert scatter_count > 1
    # Run SNP recalibrator in a scattered mode
    snp_model_path = tmp_prefix / 'snp_model'
    if not can_reuse(snp_model_path):
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
        model_j.depends_on(*jobs)
        jobs.append(model_j)
        b.write_output(model_j.model_file, str(snp_model_path))
    snp_model = b.read_input(str(snp_model_path))

    intervals_j, intervals = get_intervals(
        b=b,
        scatter_count=scatter_count,
        source_intervals_path=intervals_path,
        job_attrs=job_attrs,
        output_prefix=tmp_prefix / f'intervals_{scatter_count}',
    )
    if intervals_j:
        jobs.append(intervals_j)

    snps_recal_paths = [
        tmp_prefix / f'snp_recalibrations_{i}' for i in range(scatter_count)
    ]
    snps_tranches_paths = [
        tmp_prefix / f'snp_tranches_{i}' for i in range(scatter_count)
    ]
    scattered_jobs = []
    for idx in range(scatter_count):
        if not can_reuse(
            [
                snps_recal_paths[idx],
                str(snps_recal_paths[idx]) + '.idx',
                snps_tranches_paths[idx],
            ]
        ):
            snps_recal_j = snps_recalibrator_scattered(
                b,
                siteonly_vcf=siteonly_vcf,
                interval=intervals[idx],
                model_file=snp_model,
                hapmap_resource_vcf=resources['hapmap'],
                omni_resource_vcf=resources['omni'],
                one_thousand_genomes_resource_vcf=resources['one_thousand_genomes'],
                dbsnp_resource_vcf=resources['dbsnp'],
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                is_small_callset=is_small_callset,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            snps_recal_j.depends_on(*jobs)
            scattered_jobs.append(snps_recal_j)
            b.write_output(snps_recal_j.recalibration, str(snps_recal_paths[idx]))
            b.write_output(
                snps_recal_j.recalibration_idx, str(snps_recal_paths[idx]) + '.idx'
            )
            b.write_output(snps_recal_j.tranches, str(snps_tranches_paths[idx]))
    jobs.extend(scattered_jobs)
    snps_recalibrations = [b.read_input(str(p)) for p in snps_recal_paths]
    snps_recalibration_idxs = [b.read_input(str(p) + '.idx') for p in snps_recal_paths]
    snps_tranches = [b.read_input(str(p)) for p in snps_tranches_paths]

    snp_gathered_tranches_path = tmp_prefix / 'snp_gathered_tranches'
    if not can_reuse(snp_gathered_tranches_path):
        snps_gather_tranches_j = snps_gather_tranches_job(
            b,
            tranches=snps_tranches,
            disk_size=small_disk,
            job_attrs=job_attrs,
        )
        snps_gather_tranches_j.depends_on(*jobs)
        jobs.append(snps_gather_tranches_j)
        b.write_output(
            snps_gather_tranches_j.out_tranches, str(snp_gathered_tranches_path)
        )
    snp_gathered_tranches = b.read_input(str(snp_gathered_tranches_path))

    interval_snps_applied_vcf_paths = [
        tmp_prefix / f'interval_snps_applied_{idx}.vcf.gz'
        for idx in range(scatter_count)
    ]
    scattered_apply_jobs = []
    for idx in range(scatter_count):
        if not can_reuse(
            [
                interval_snps_applied_vcf_paths[idx],
                to_path(str(interval_snps_applied_vcf_paths[idx]) + '.tbi'),
            ]
        ):
            j = apply_recalibration_snps(
                b,
                input_vcf=siteonly_vcf,
                interval=intervals[idx],
                recalibration=snps_recalibrations[idx],
                recalibration_idx=snps_recalibration_idxs[idx],
                tranches=snp_gathered_tranches,
                disk_size=huge_disk,
                use_as_annotations=use_as_annotations,
                filter_level=snp_filter_level,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            j.depends_on(*jobs)
            scattered_apply_jobs.append(j)
            b.write_output(j.output_vcf, str(interval_snps_applied_vcf_paths[idx]))
            b.write_output(
                j.output_tbi, str(interval_snps_applied_vcf_paths[idx]) + '.tbi'
            )
    jobs.extend(scattered_apply_jobs)
    interval_snps_applied_vcfs = [
        b.read_input_group(
            **{
                'vcf.gz': str(p),
                'vcf.gz.tbi': str(p) + '.tbi',
            }
        )
        for p in interval_snps_applied_vcf_paths
    ]

    gathered_vcf_path = tmp_prefix / 'gathered.vcf.gz'
    snps_applied_gathered_jobs = gather_vcfs(
        b=b,
        input_vcfs=interval_snps_applied_vcfs,
        site_only=True,
        sample_count=gvcf_count,
        job_attrs=job_attrs,
        out_vcf_path=gathered_vcf_path,
    )
    for j in snps_applied_gathered_jobs:
        j.name = f'VQSR: {j.name}'
        j.depends_on(*jobs)
        jobs.append(j)
    snps_applied_gathered_vcf = b.read_input_group(
        **{
            'vcf.gz': str(gathered_vcf_path),
            'vcf.gz.tbi': f'{gathered_vcf_path}.tbi',
        }
    )

    apply_indel_j = apply_recalibration_indels(
        b,
        input_vcf=snps_applied_gathered_vcf,
        recalibration=indel_recalibration,
        recalibration_idx=indel_recalibration_idx,
        tranches=indel_tranches,
        disk_size=small_disk,
        use_as_annotations=use_as_annotations,
        filter_level=indel_filter_level,
        job_attrs=job_attrs,
    )
    apply_indel_j.depends_on(*jobs)
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

    gatk --java-options \
      "-Xms{res.get_java_mem_mb()}m \
      -XX:+UseParallelGC \
      -XX:ParallelGCThreads={res.get_nthreads() - 2}" \\
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
    
    mv {j.recalibration}.idx {j.recalibration_idx}
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

    gatk --java-options \
      "-Xms{res.get_java_mem_mb()}m \
      -XX:+UseParallelGC \
      -XX:ParallelGCThreads={res.get_nthreads() - 2}" \\
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
        res = HIGHMEM.set_resources(j, ncpu=4, storage_gb=disk_size)

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

    gatk --java-options \
      "-Xms{res.get_java_mem_mb()}m \
      -XX:+UseParallelGC \
      -XX:ParallelGCThreads={res.get_nthreads() - 2}" \\
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
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 {dbsnp_resource_vcf.base}
    
    mv {j.recalibration}.idx {j.recalibration_idx}  
    """
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

    gatk --java-options \
      "-Xms{res.get_java_mem_mb()}m \
      -XX:+UseParallelGC \
      -XX:ParallelGCThreads={res.get_nthreads() - 2}" \\
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
    tranches: Sequence[hb.ResourceFile],
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
    recalibration: hb.Resource,
    recalibration_idx: hb.Resource,
    tranches: hb.ResourceFile,
    disk_size: int,
    use_as_annotations: bool,
    filter_level: float,
    interval: hb.Resource | None = None,
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

    cmd = f"""
    cp {recalibration} $BATCH_TMPDIR/recalibration
    cp {recalibration_idx} $BATCH_TMPDIR/recalibration.idx
    
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    ApplyVQSR \\
    -O $BATCH_TMPDIR/output.vcf.gz \\
    -V {input_vcf['vcf.gz']} \\
    --recal-file $BATCH_TMPDIR/recalibration \\
    --tranches-file {tranches} \\
    --truth-sensitivity-filter-level {filter_level} \\
    {f'-L {interval} ' if interval else ''} \\
    {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
    -mode SNP

    tabix -p vcf -f $BATCH_TMPDIR/output.vcf.gz
    mv $BATCH_TMPDIR/output.vcf.gz {j.output_vcf}
    mv $BATCH_TMPDIR/output.vcf.gz.tbi {j.output_tbi}
    """
    j.command(command(cmd, monitor_space=True))
    return j


def apply_recalibration_indels(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    recalibration: hb.Resource,
    recalibration_idx: hb.Resource,
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
    mv {recalibration} $BATCH_TMPDIR/recalibration
    mv {recalibration_idx} $BATCH_TMPDIR/recalibration.idx
    
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    ApplyVQSR \\
    --tmp-dir $BATCH_TMPDIR \\
    -O {j.output_vcf['vcf.gz']} \\
    -V {input_vcf['vcf.gz']} \\
    --recal-file $BATCH_TMPDIR/recalibration \\
    --tranches-file {tranches} \\
    --truth-sensitivity-filter-level {filter_level} \\
    {f'-L {interval} ' if interval else ''} \\
    {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
    -mode INDEL

    tabix -p vcf -f {j.output_vcf['vcf.gz']}
    """
    j.command(command(cmd, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j
