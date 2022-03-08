"""
Create Hail Batch jobs for variant calling in individual samples.
"""

import logging
from os.path import join, basename
from typing import Optional, List, Tuple

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import images, ref_data, buckets
from cpg_pipes.jobs import split_intervals
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.pipeline.analysis import CramPath
from cpg_pipes.pipeline.smdb import SMDB

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def produce_gvcf(
    b: hb.Batch,
    sample_name: str,
    dataset_name: str,
    tmp_bucket: str,
    cram: CramPath,
    output_path: Optional[str] = None,
    number_of_intervals: int = 1,
    intervals: Optional[hb.ResourceGroup] = None,
    overwrite: bool = True,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    dragen_mode: bool = False,
) -> Job:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    if output_path and buckets.can_reuse(output_path, overwrite):
        return b.new_job('Make GVCF [reuse]', dict(sample=sample_name, dataset=dataset_name))

    depends_on = depends_on or []

    first_j, last_j, hc_gvcf_path = haplotype_caller(
        b=b,
        sample_name=sample_name,
        dataset_name=dataset_name,
        tmp_bucket=tmp_bucket,
        cram=cram,
        number_of_intervals=number_of_intervals, 
        intervals=intervals,
        overwrite=overwrite, 
        depends_on=depends_on,
        dragen_mode=dragen_mode,
    )

    postproc_j = postproc_gvcf(
        b=b,
        sample_name=sample_name,
        dataset_name=dataset_name,
        gvcf_path=hc_gvcf_path,
        output_path=output_path,
        overwrite=overwrite,
        depends_on=depends_on + [last_j],
    )
    last_j = postproc_j

    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type='gvcf',
            output_path=output_path,
            sample_names=[sample_name],
            dataset_name=dataset_name,
            first_j=first_j,
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def haplotype_caller(
    b: hb.Batch,
    sample_name: str,
    dataset_name: str,
    tmp_bucket: str,
    cram: CramPath,
    number_of_intervals: int = 1,
    intervals: Optional[hb.ResourceGroup] = None,
    overwrite: bool = True,
    depends_on: Optional[List[Job]] = None,
    dragen_mode: bool = False,
) -> Tuple[Job, Job, str]:
    """
    Run haplotype caller in parallel sharded by intervals. 
    Returns the first and the last job object, and path to the output GVCF file.
    """
    hc_gvcf_path = join(tmp_bucket, 'haplotypecaller', f'{sample_name}.g.vcf.gz')
    if buckets.can_reuse(hc_gvcf_path, overwrite):
        first_j = last_j = b.new_job('HaplotypeCaller [reuse]', dict(
            sample=sample_name, dataset=dataset_name))
        return first_j, last_j, hc_gvcf_path

    hc_jobs = []
    if number_of_intervals > 1:
        if intervals is None:
            intervals = split_intervals.get_intervals(
                b=b,
                scatter_count=number_of_intervals,
                out_bucket=join(tmp_bucket, 'intervals'),
            )

        # Splitting variant calling by intervals
        for idx in range(number_of_intervals):
            hc_jobs.append(
                _haplotype_caller_one(
                    b,
                    sample_name=sample_name,
                    dataset_name=dataset_name,
                    cram=cram,
                    interval=intervals[f'interval_{idx}'],
                    interval_idx=idx,
                    number_of_intervals=number_of_intervals,
                    depends_on=depends_on,
                    dragen_mode=dragen_mode,
                    overwrite=overwrite,
                )
            )
        last_j = merge_gvcfs_job(
            b=b,
            sample_name=sample_name,
            dataset_name=dataset_name,
            gvcfs=[j.output_gvcf for j in hc_jobs],
            out_gvcf_path=hc_gvcf_path,
            overwrite=overwrite,
        )
    else:
        hc_j = _haplotype_caller_one(
            b,
            sample_name=sample_name,
            dataset_name=dataset_name,
            cram=cram,
            depends_on=depends_on,
            out_gvcf_path=hc_gvcf_path,
            overwrite=overwrite,
        )
        hc_jobs.append(hc_j)
        last_j = hc_j

    first_j = hc_jobs[0]
    if depends_on:
        first_j.depends_on(*depends_on)
    
    return first_j, last_j, hc_gvcf_path


def _haplotype_caller_one(
    b: hb.Batch,
    sample_name: str,
    dataset_name: str,
    cram: CramPath,
    interval: Optional[hb.ResourceFile] = None,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
    out_gvcf_path: Optional[str] = None,
    overwrite: bool = True,
    dragen_mode: bool = False,
) -> Job:
    """
    Add one HaplotypeCaller job on an interval
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name, dict(sample=sample_name, dataset=dataset_name))
    if buckets.can_reuse(out_gvcf_path, overwrite):
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
    if depends_on:
        j.depends_on(*depends_on)

    ref_fasta = ref_data.REF_FASTA
    ref_fai = ref_data.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )
    
    # cram = b.read_input_group(**{
    #     'cram': cram_fpath,
    #     'cram.crai': crai_fpath or (cram_fpath + '.crai'),
    # })

    cmd = f"""\
    CRAM=/io/batch/{sample_name}.cram
    CRAI=/io/batch/{sample_name}.cram.crai

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {cram.path} $CRAM
    retry_gs_cp {cram.index_path} $CRAI

    # Copying reference data as well to avoid crazy logging costs 
    # for region requests
    retry_gs_cp {ref_fasta} /io/batch/{basename(ref_fasta)}
    retry_gs_cp {ref_fai}   /io/batch/{basename(ref_fai)}
    retry_gs_cp {ref_dict}  /io/batch/{basename(ref_dict)}

    gatk --java-options "-Xms{job_res.get_java_mem_mb()}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R /io/batch/{basename(ref_fasta)} \\
    -I $CRAM \\
    --read-index $CRAI \\
    {f"-L {interval} " if interval is not None else ""} \\
    --disable-spanning-event-genotyping \\
    {"--dragen-mode " if dragen_mode else ""} \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 20 \\
    -ERC GVCF
    """
    j.command(wrap_command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    dataset_name: str,
    gvcfs: List[hb.ResourceGroup],
    out_gvcf_path: Optional[str],
    overwrite: bool = True,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """
    job_name = f'Merge {len(gvcfs)} GVCFs'
    j = b.new_job(job_name, dict(sample=sample_name, dataset=dataset_name))
    if buckets.can_reuse(out_gvcf_path, overwrite):
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
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def postproc_gvcf(
    b: hb.Batch,
    gvcf_path: str,
    sample_name: str,
    dataset_name: str,
    overwrite: bool,
    output_path: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    1. Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration, and reduce the number of GVCF blocking bins to 2.
    2. Subsets GVCF to main, not-alt chromosomes to avoid downstream errors.
    3. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    4. Renames the GVCF sample name to use CPG ID.
    """
    if output_path and buckets.can_reuse(output_path, overwrite):
        return b.new_job('Postproc GVCF [reuse]', dict(sample=sample_name, dataset=dataset_name))

    logger.info(f'Adding GVCF postproc job for sample {sample_name}, gvcf {gvcf_path}')

    j = b.new_job(f'ReblockGVCF', dict(sample=sample_name, dataset=dataset_name))
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

    ref_fasta = ref_data.REF_FASTA
    ref_fai = ref_data.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )

    cmd = f"""\
    GVCF=/io/batch/{sample_name}.g.vcf.gz
    GVCF_NODP=/io/batch/{sample_name}-nodp.g.vcf.gz
    REBLOCKED=/io/batch/{sample_name}-reblocked.g.vcf.gz

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {gvcf_path} $GVCF
    retry_gs_cp {ref_data.NOALT_REGIONS} /io/batch/noalt-regions.bed

    # Copying reference data as well to avoid crazy logging costs 
    # for region requests
    retry_gs_cp {ref_fasta} /io/batch/{basename(ref_fasta)}
    retry_gs_cp {ref_fai}   /io/batch/{basename(ref_fai)}
    retry_gs_cp {ref_dict}  /io/batch/{basename(ref_dict)}

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

    gatk --java-options "-Xms{job_res.get_java_mem_mb()}g" \\
    ReblockGVCF \\
    --reference /io/batch/{basename(ref_fasta)} \\
    -V $GVCF_NODP \\
    -do-qual-approx \\
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools query -l $GVCF)

    bcftools view $REBLOCKED -T /io/batch/noalt-regions.bed \\
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
