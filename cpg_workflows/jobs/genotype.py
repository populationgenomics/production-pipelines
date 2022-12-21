"""
CRAM to GVCF: create Hail Batch jobs to genotype individual samples.
"""

import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path, command

from cpg_workflows.filetypes import CramPath, GvcfPath
from cpg_workflows.resources import STANDARD, HIGHMEM
from cpg_workflows import utils
from .picard import get_intervals


def genotype(
    b: hb.Batch,
    sample_name: str,
    tmp_prefix: Path,
    cram_path: CramPath,
    job_attrs: dict[str, str] | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> list[Job]:
    """
    Takes a CRAM file and runs GATK tools to make a GVCF.
    """
    hc_gvcf_path = tmp_prefix / 'haplotypecaller' / f'{sample_name}.g.vcf.gz'
    if utils.can_reuse(output_path, overwrite):
        return []

    jobs = haplotype_caller(
        b=b,
        sample_name=sample_name,
        job_attrs=job_attrs,
        output_path=hc_gvcf_path,
        cram_path=cram_path,
        tmp_prefix=tmp_prefix,
        scatter_count=50,
        overwrite=overwrite,
        dragen_mode=dragen_mode,
    )
    postproc_j = postproc_gvcf(
        b=b,
        gvcf_path=GvcfPath(hc_gvcf_path),
        sample_name=sample_name,
        job_attrs=job_attrs,
        output_path=output_path,
        overwrite=overwrite,
    )
    postproc_j.depends_on(*jobs)
    return jobs + [postproc_j]


intervals: list[hb.ResourceFile] | None = None


def haplotype_caller(
    b: hb.Batch,
    sample_name: str,
    cram_path: CramPath,
    tmp_prefix: Path,
    scatter_count: int,
    job_attrs: dict[str, str] | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> list[Job]:
    """
    Run GATK Haplotype Caller in parallel, split by intervals.
    """
    if utils.can_reuse(output_path, overwrite):
        return [b.new_job('HaplotypeCaller [reuse]', job_attrs)]

    jobs: list[Job] = []

    if scatter_count > 1:
        global intervals
        if intervals is None:
            intervals_j, intervals = get_intervals(
                b=b,
                scatter_count=scatter_count,
                source_intervals_path=get_config()['workflow'].get('intervals_path'),
                job_attrs=job_attrs,
                output_prefix=tmp_prefix / f'intervals_{scatter_count}',
            )
            if intervals_j:
                jobs.append(intervals_j)

        hc_jobs = []
        # Splitting variant calling by intervals
        for idx in range(scatter_count):
            assert intervals[idx], intervals
            j = _haplotype_caller_one(
                b,
                sample_name=sample_name,
                cram_path=cram_path,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
                interval=intervals[idx],
                dragen_mode=dragen_mode,
                overwrite=overwrite,
            )
            hc_jobs.append(j)
        jobs.extend(hc_jobs)
        merge_j = merge_gvcfs_job(
            b=b,
            sample_name=sample_name,
            gvcf_groups=[j.output_gvcf for j in hc_jobs],
            job_attrs=job_attrs,
            out_gvcf_path=output_path,
            overwrite=overwrite,
        )
        jobs.append(merge_j)
    else:
        hc_j = _haplotype_caller_one(
            b,
            sample_name=sample_name,
            job_attrs=job_attrs,
            cram_path=cram_path,
            out_gvcf_path=output_path,
            overwrite=overwrite,
        )
        jobs.append(hc_j)

    return jobs


def _haplotype_caller_one(
    b: hb.Batch,
    sample_name: str,
    cram_path: CramPath,
    job_attrs: dict | None = None,
    interval: hb.Resource | None = None,
    out_gvcf_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> Job:
    """
    Add one GATK HaplotypeCaller job on an interval.
    """
    job_name = 'HaplotypeCaller'
    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='gatk HaplotypeCaller'))
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name = f'{j.name} [reuse]'
        return j

    j.image(image_path('gatk'))

    # Enough storage to localize CRAMs (can't pass GCS URL to CRAM to gatk directly
    # because we will hit GCP egress bandwidth limit:
    # https://batch.hail.populationgenomics.org.au/batches/7493/jobs/2)
    # CRAMs can be as big as 80G:
    # https://batch.hail.populationgenomics.org.au/batches/74042/jobs/3346
    # plus we need enough space to fit output GVCF (1G) and reference data (5G).
    # HaplotypeCaller is not parallelized, so we request the minimal possible chunk
    # of a Hail worker (2 cores), plus with the `highmem` instance that would
    # give enough memory: 13G. That's not going to give enough disk storage, so we
    # are explicitly requesting more storage.
    storage_gb = None  # avoid extra disk by default
    if get_config()['workflow']['sequencing_type'] == 'genome':
        storage_gb = 100
    job_res = HIGHMEM.request_resources(ncpu=2)
    # enough for input CRAM and output GVCF
    job_res.attach_disk_storage_gb = storage_gb
    job_res.set_to_job(j)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    reference = fasta_res_group(b)

    assert isinstance(j.output_gvcf, hb.ResourceGroup)

    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{sample_name}.cram
    CRAI=$BATCH_TMPDIR/{sample_name}.cram.crai

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    gatk --java-options \
    "-Xms{job_res.get_java_mem_mb()}m \
    -Xmx{job_res.get_java_mem_mb()}m \
    -XX:GCTimeLimit=50 \
    -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R {reference.base} \\
    -I $CRAM \\
    --read-index $CRAI \\
    {f"-L {interval} " if interval is not None else ""} \\
    --disable-spanning-event-genotyping \\
    {"--dragen-mode " if dragen_mode else ""} \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
    -ERC GVCF \\
    --create-output-variant-index
    """
    j.command(
        command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
    )
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


def merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    gvcf_groups: list[hb.Resource],
    job_attrs: dict | None = None,
    out_gvcf_path: Path | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample-wide GVCF file.
    """
    job_name = f'Merge {len(gvcf_groups)} GVCFs'
    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='picard MergeVcfs'))
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name = f'{j.name} [reuse]'
        return j

    j.image(image_path('picard'))
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcf_groups) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ''
    for gvcf_group in gvcf_groups:
        assert isinstance(gvcf_group, hb.ResourceGroup)
        input_cmd += f'INPUT={gvcf_group["g.vcf.gz"]} '

    assert isinstance(j.output_gvcf, hb.ResourceGroup)
    cmd = f"""
    picard -Xms{java_mem}g \
    MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}
    """
    j.command(command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


def postproc_gvcf(
    b: hb.Batch,
    gvcf_path: GvcfPath,
    sample_name: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
    output_path: Path | None = None,
    depends_on: list[Job] | None = None,
) -> Job:
    """
    1. Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration, and reduce the number of GVCF blocking bins to 2.
    2. Subsets GVCF to main, not-alt chromosomes to avoid downstream errors.
    3. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    4. Renames the GVCF sample name to use CPG ID.
    """
    logging.info(f'Adding GVCF postproc job for sample {sample_name}, gvcf {gvcf_path}')
    job_name = 'Postproc GVCF'
    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='gatk ReblockGVCF'))
    if utils.can_reuse(output_path, overwrite):
        return b.new_job(job_name + ' [reuse]', job_attrs)

    j.image(image_path('gatk'))

    # We need at least 2 CPU, so on 16-core instance it would be 8 jobs,
    # meaning we have more than enough disk (265/8=33.125G).
    # Enough to fit a pre-reblocked GVCF, which can be as big as 10G,
    # the reblocked result (1G), and ref data (5G).
    job_res = STANDARD.set_resources(j, ncpu=2, storage_gb=20)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    reference = fasta_res_group(b)
    noalt_regions = b.read_input(str(reference_path('broad/noalt_bed')))
    gvcf = b.read_input(str(gvcf_path.path))
    gq_bands = get_config()['workflow']['reblock_gq_bands']

    assert isinstance(j.output_gvcf, hb.ResourceGroup)

    cmd = f"""\
    GVCF={gvcf}
    GVCF_NODP=$BATCH_TMPDIR/{sample_name}-nodp.g.vcf.gz
    REBLOCKED=$BATCH_TMPDIR/{sample_name}-reblocked.g.vcf.gz

    # Reindexing just to make sure the index is not corrupted
    bcftools index --tbi $GVCF

    # Remove INFO/DP field, which contradicts the FORMAT/DP, in the way that
    # it has _all_ reads, not just variant-calling-usable reads. If we keep INFO/DP,
    # ReblockGVCF would prioritize it over FORMAT/DP and change FORMAT/DP to INFO/DP
    # in the resulting merged blocks. It would pick the highest INFO/DP when merging
    # multiple blocks, so a variant in a small homopolymer region (surrounded by
    # long DP=0 areas), that attracted piles of low-MQ reads with INFO/DP=1000
    # will translate into a long GQ<20 block with the same FORMAT/DP=1000, 
    # which is wrong, because most of this block has no reads.
    bcftools view $GVCF \\
    | bcftools annotate -x INFO/DP \\
    | bcftools view -Oz -o $GVCF_NODP
    tabix -p vcf $GVCF_NODP

    gatk --java-options "-Xms{job_res.get_java_mem_mb()}m" \\
    ReblockGVCF \\
    --reference {reference.base} \\
    -V $GVCF_NODP \\
    -do-qual-approx \\
    --floor-blocks {' '.join(f'-GQB {b}' for b in gq_bands)} \\
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools query -l $GVCF)

    bcftools view $REBLOCKED -T {noalt_regions} \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "$EXISTING_SN {sample_name}") \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    tabix -p vcf {j.output_gvcf['g.vcf.gz']}
    """
    j.command(
        command(cmd, setup_gcp=True, monitor_space=True, define_retry_function=True)
    )
    if output_path:
        b.write_output(j.output_gvcf, str(output_path).replace('.g.vcf.gz', ''))
    if depends_on:
        j.depends_on(*depends_on)
    return j
