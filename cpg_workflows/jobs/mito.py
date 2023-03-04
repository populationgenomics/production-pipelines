"""
Create Hail Batch jobs to run STRipy
"""
import hailtop.batch as hb
from hailtop.batch.job import Job
from textwrap import dedent

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse
from cpg_workflows.targets import Sample

from cpg_workflows.jobs import picard


def subset_cram_to_chrM(
    b,
    cram_path: CramPath,
    output_bam_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Extract read pairs with at least one read mapping to chrM

    This uses the same extraction command as the broad pipeline
    as I am not clear the specifics *may* be important to how
    NUMT mapping reads are handled
    """

    # If we can resuse existing output, return the existing cram
    if can_reuse(
        [
            output_bam_path,
        ],
        overwrite,
    ):
        return None, b.read_input_group(
            **{
                'bam': str(output_bam_path),
                'bam.bai': str(output_bam_path) + '.bai',
            }
        )

    job_attrs = job_attrs or {}
    j = b.new_job('subset_cram_to_chrM', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    reference = fasta_res_group(b)
    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bai',
        }
    )

    cmd = f"""
        CRAM=$BATCH_TMPDIR/{cram_path.path.name}
        CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path.path)} $CRAM
        retry_gs_cp {str(cram_path.index_path)} $CRAI

        gatk PrintReads \
            -R {reference.base} \
            -L chrM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -I $CRAM \
            --read-index $CRAI \
            -O {j.output_bam.bam}

        dirname {j.output_bam.bam}
        ls -l $(dirname {j.output_bam.bam})
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_bam, str(output_bam_path.with_suffix('')))

    return j, j.output_bam


def mito_realign(
    b,
    sample_id: str,
    input_bam: hb.ResourceGroup,
    output_cram_path: Path,
    mito_ref: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Re-align reads to mito genome
    """
    if can_reuse(
        [
            output_cram_path,
        ],
        overwrite,
    ):
        return []

    job_attrs = job_attrs or {}
    j = b.new_job('mito_realign', job_attrs)
    j.image(image_path('bwa'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    nthreads = res.get_nthreads()

    cmd = f"""\
        bazam -Xmx16g \
            -n{min(nthreads, 6)} -bam {input_bam.bam} | \
        bwa \
            mem -K 100000000 -p -v 3 -t 2 -Y {mito_ref.fasta} \
            -R '@RG\\tID:{sample_id}\\tSM:{sample_id}' \
            - | \
        samtools view -bSu - | \
        samtools sort -o {j.raw_cram}
        """
    j.command(command(cmd, define_retry_function=True))

    mkdup_j = picard.markdup(
        b=b,
        sorted_bam=j.raw_cram,
        output_path=output_cram_path,
        out_markdup_metrics_path=output_cram_path.with_suffix('.markduplicates-metrics')
    )

    return [j, mkdup_j]


def mito_mutect2(
    b,
    cram_path: Path,
    reference: hb.ResourceGroup,
    region: str,
    max_reads_per_alignment_start: int = 75,
    job_attrs: dict | None = None,
) -> Job:
    """
    Call mutect2 on a mito genome

    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_mutect2', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    java_mem_mb = res.get_java_mem_mb()

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        CRAM=$BATCH_TMPDIR/{cram_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path)} $CRAM
        retry_gs_cp {str(cram_path)}.crai $CRAM.crai

        # We need to create these files regardless, even if they stay empty
        # touch bamout.bam

        gatk --java-options "-Xmx{java_mem_mb}m" Mutect2 \
            -R {reference.fasta} \
            -I $CRAM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {j.output_vcf['vcf.gz']} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            -L {region}
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def liftover_and_combine_vcfs(
    b,
    vcf: hb.ResourceGroup,
    shifted_vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    shift_back_chain: hb.ResourceFile,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Lifts over shifted vcf of control region and combines it with the rest of the chrM calls.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('liftover_and_combine_vcfs', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        picard LiftoverVcf \
        I={shifted_vcf['vcf.gz']} \
        O={j.lifted_vcf} \
        R={reference.fasta} \
        CHAIN={shift_back_chain} \
        REJECT={j.rejected_vcf}

        picard MergeVcfs \
        I={vcf['vcf.gz']} \
        I={j.lifted_vcf} \
        O={j.rejected_vcf}.merged.vcf
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def filter_variants(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Mutect2 Filtering for calling Snps and Indels
    """
    job_attrs = job_attrs or {}
    j = b.new_job('filter_variants', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
      gatk --java-options "-Xmx2500m" FilterMutectCalls -V {vcf['vcf.gz']} \
        -R {reference.fasta} \
        -O filtered.vcf \
        --stats {raw_vcf_stats} \
        {m2_extra_filtering_args} \
        --max-alt-allele-count {max_alt_allele_count} \
        --mitochondria-mode \
        --min-allele-fraction  {vaf_filter_threshold} \
        --f-score-beta  {f_score_beta} \
        --contamination-estimate  {max_contamination}

      gatk VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        --mask {blacklisted_sites} \
        --mask-name "blacklisted_site"

    """

    j.command(command(cmd, define_retry_function=True))

    return j


def genotype_mito(
    b,
    cram_path: Path,
    shifted_cram_path: Path,
    output_vcf_path: Path,
    mito_reff: hb.ResourceGroup,
    shifted_mito_reff: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job | None]:
    """
    Genotype mito genome using mutect2

    """

    # Call variants
    call_j = mito_mutect2(
        b=b,
        cram_path=cram_path,
        reference=mito_reff,
        region='chrM:576-16024',  # Exclude the control region.
        job_attrs=job_attrs
    )

    shifted_call_j = mito_mutect2(
        b=b,
        cram_path=shifted_cram_path,
        reference=shifted_mito_reff,
        region='chrM:8025-9144',  # Only call inside the control region.
        job_attrs=job_attrs

    )
    jobs = [call_j, shifted_call_j]

    # Merge two VCFs
    merge_j = liftover_and_combine_vcfs(
        b=b,
        vcf=call_j.output_vcf,
        shifted_vcf=shifted_call_j.output_vcf,
        reference=mito_reff,  # Does this need to be full genome ref?
        shift_back_chain=shifted_mito_reff.shift_back_chain,
        job_attrs=job_attrs
    )
    jobs.append(merge_j)

    return jobs
