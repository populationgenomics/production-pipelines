"""
Create Hail Batch jobs for variant calling in individual samples.
"""

import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes.types import CramPath, GvcfPath, SequencingType
from cpg_pipes import images, utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs import split_intervals
from cpg_pipes.refdata import RefData

logger = logging.getLogger(__file__)


def produce_gvcf(
    b: hb.Batch,
    sample_name: str,
    sequencing_type: SequencingType,
    tmp_bucket: Path,
    cram_path: CramPath,
    refs: RefData,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    number_of_intervals: int = 1,
    intervals: list[hb.Resource] | None = None,
    overwrite: bool = True,
    dragen_mode: bool = False,
) -> list[Job]:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    if utils.can_reuse(output_path, overwrite):
        return [b.new_job('Make GVCF [reuse]', job_attrs)]

    hc_gvcf_path = tmp_bucket / 'haplotypecaller' / f'{sample_name}.g.vcf.gz'

    jobs = haplotype_caller(
        b=b,
        sample_name=sample_name,
        sequencing_type=sequencing_type,
        refs=refs,
        job_attrs=job_attrs,
        output_path=hc_gvcf_path,
        tmp_bucket=tmp_bucket,
        cram_path=cram_path,
        number_of_intervals=number_of_intervals,
        intervals=intervals,
        overwrite=overwrite,
        dragen_mode=dragen_mode,
    )

    postproc_j = postproc_gvcf(
        b=b,
        gvcf_path=GvcfPath(hc_gvcf_path),
        sample_name=sample_name,
        refs=refs,
        job_attrs=job_attrs,
        output_path=output_path,
        overwrite=overwrite,
    )
    postproc_j.depends_on(*jobs)
    
    return jobs + [postproc_j]


def haplotype_caller(
    b: hb.Batch,
    sample_name: str,
    sequencing_type: SequencingType,
    tmp_bucket: Path,
    cram_path: CramPath,
    refs: RefData,
    job_attrs: dict | None = None,
    output_path: Path | None = None,
    number_of_intervals: int = 1,
    intervals: list[hb.Resource] | None = None,
    overwrite: bool = True,
    dragen_mode: bool = False,
) -> list[Job]:
    """
    Run haplotype caller in parallel sharded by intervals. 
    Returns jobs and path to the output GVCF file.
    """
    if utils.can_reuse(output_path, overwrite):
        return [b.new_job('HaplotypeCaller [reuse]', job_attrs)]
    
    jobs: list[Job] = []
    if number_of_intervals > 1:
        if intervals is None:
            intervals_j = split_intervals.get_intervals(
                b=b,
                refs=refs,
                sequencing_type=sequencing_type,
                scatter_count=number_of_intervals,
                out_bucket=tmp_bucket / 'intervals',
            )
            jobs.append(intervals_j)
            intervals = [intervals_j[f'intervals{i}.list'] for i in range(number_of_intervals)]

        hc_jobs = []
        # Splitting variant calling by intervals
        for idx in range(number_of_intervals):
            j = _haplotype_caller_one(
                b,
                sample_name=sample_name,
                cram_path=cram_path,
                refs=refs,
                job_attrs=job_attrs,
                interval=intervals[idx],
                interval_idx=idx,
                number_of_intervals=number_of_intervals,
                dragen_mode=dragen_mode,
                overwrite=overwrite,
            )
            hc_jobs.append(j)
        merge_j = merge_gvcfs_job(
            b=b,
            sample_name=sample_name,
            job_attrs=job_attrs,
            gvcfs=[j.output_gvcf for j in hc_jobs],
            out_gvcf_path=output_path,
            overwrite=overwrite,
        )
        jobs.extend(jobs + [merge_j])
    else:
        hc_j = _haplotype_caller_one(
            b,
            sample_name=sample_name,
            refs=refs,
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
    refs: RefData,
    job_attrs: dict | None = None,
    interval: hb.Resource | None = None,
    interval_idx: int | None = None,
    number_of_intervals: int = 1,
    out_gvcf_path: Path | None = None,
    overwrite: bool = True,
    dragen_mode: bool = False,
) -> Job:
    """
    Add one HaplotypeCaller job on an interval
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name, job_attrs)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.GATK_IMAGE)

    # Enough storage to localize CRAMs (can't pass GCS URL to CRAM to gatk directly
    # because we will hit GCP egress bandwidth limit:
    # https://batch.hail.populationgenomics.org.au/batches/7493/jobs/2)
    # 45 should be enough to fit a CRAM (30G), GVCF (1G), and ref data (5G),
    # and we can squeeze 4 jobs on a 32-core machine (185G/4=46.25G) or
    # 5 jobs on a 16-core machine (265G/4=53G)
    job_res = STANDARD.set_resources(j, storage_gb=45)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    reference = refs.fasta_res_group(b)
    
    cmd = f"""\
    CRAM=/io/batch/{sample_name}.cram
    CRAI=/io/batch/{sample_name}.cram.crai

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
    -GQB 20 \\
    -ERC GVCF \\
    --create-output-variant-index
    """
    j.command(wrap_command(
        cmd, monitor_space=True, setup_gcp=True, define_retry_function=True
    ))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


def merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    gvcfs: list[hb.ResourceGroup],
    job_attrs: dict | None = None,
    out_gvcf_path: Path | None = None,
    overwrite: bool = True,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """
    job_name = f'Merge {len(gvcfs)} GVCFs'
    j = b.new_job(job_name, job_attrs)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j
    
    j.image(images.SAMTOOLS_PICARD_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)
    cmd = f"""
    picard -Xms{java_mem}g \
    MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


def postproc_gvcf(
    b: hb.Batch,
    gvcf_path: GvcfPath,
    sample_name: str,
    refs: RefData,
    job_attrs: dict | None = None,
    overwrite: bool = True,
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
    if utils.can_reuse(output_path, overwrite):
        return b.new_job('Postproc GVCF [reuse]', job_attrs)

    logger.info(f'Adding GVCF postproc job for sample {sample_name}, gvcf {gvcf_path}')

    j = b.new_job(f'ReblockGVCF', job_attrs)
    j.image(images.GATK_IMAGE)

    # Enough to fit a pre-reblocked GVCF, which can be as big as 10G,
    # the reblocked result (1G), and ref data (5G). 
    # We need at least 2 CPU, so on 16-core instance it would be 8 jobs,
    # meaning we have more than enough disk (265/8=33.125G). 
    job_res = STANDARD.set_resources(j, storage_gb=20)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    reference = refs.fasta_res_group(b)
    noalt_regions = b.read_input(str(refs.noalt_regions))
    gvcf = b.read_input(str(gvcf_path.path))

    cmd = f"""\
    GVCF={gvcf}
    GVCF_NODP=/io/batch/{sample_name}-nodp.g.vcf.gz
    REBLOCKED=/io/batch/{sample_name}-reblocked.g.vcf.gz

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
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools query -l $GVCF)

    bcftools view $REBLOCKED -T {noalt_regions} \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "$EXISTING_SN {sample_name}") \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    tabix -p vcf {j.output_gvcf['g.vcf.gz']}
    """
    j.command(wrap_command(
        cmd, setup_gcp=True, monitor_space=True, define_retry_function=True
    ))
    if output_path:
        b.write_output(j.output_gvcf, str(output_path).replace('.g.vcf.gz', ''))
    if depends_on:
        j.depends_on(*depends_on)
    return j
