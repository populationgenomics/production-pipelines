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
            'bam.bai': '{root}.bam.bai',
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
        ls -l $(dirname {j.output_bam.cram})
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_bam, str(output_bam_path.with_suffix('')))

    return j, j.output_bam


def mito_realign(
    b,
    input_cram: hb.ResourceGroup,
    output_cram_path: Path,
    shifted: bool,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Re-align reads to mito genome
    """
    # if can_reuse(
    #     [
    #         output_cram_path,
    #     ],
    #     overwrite,
    # ):
    #     return None

    job_attrs = job_attrs or {}
    j = b.new_job('mito_realign', job_attrs)
    j.image(image_path('bwa'))

    reference = fasta_res_group(b)
    if shifted:
        mt_ref = b.read_input_group(
            dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict',
            fasta='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta',
            amb='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb',
            ann='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann',
            bwt='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt',
            fai='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai',
            pac='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac',
            sa='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa',
        )
    else:
        mt_ref = b.read_input_group(
            dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict',
            fasta='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta',
            amb='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb',
            ann='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann',
            bwt='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt',
            fai='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai',
            pac='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac',
            sa='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa',
        )

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    nthreads = res.get_nthreads()

    cmd = dedent(
        f"""\
        bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
            -n{min(nthreads, 6)} -bam {input_cram.cram} -L chrM | \
         bwa mem -K 100000000 -p -v 3 -t 2 -Y {mt_ref.fasta} - | \
         samtools view -bSu - | \
        samtools sort -o {j.raw_cram}
        """
    )

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
    cram_path: CramPath,
    reference: hb.ResourceGroup,
    max_reads_per_alignment_start: int = 75,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Call mutect2 on a mito genome

    """
    job_attrs = job_attrs or {}
    j = b.new_job('genotype_mito', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    java_mem_mb = res.get_java_mem_mb()

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        CRAM=$BATCH_TMPDIR/{cram_path.path.name}
        CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path.path)} $CRAM
        retry_gs_cp {str(cram_path.index_path)} $CRAI

        # We need to create these files regardless, even if they stay empty
        # touch bamout.bam

        gatk --java-options "-Xmx~{java_mem_mb}m" Mutect2 \
            -R {reference.fasta} \
            -I $CRAM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O ~{output_vcf.vcf.gz} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {max_reads_per_alignment_start} \
            --max-mnp-distance 0
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_cram, str(output_cram_path.with_suffix('')))

    return j


def liftover_and_combine_vcfs(
    b,
    vcf: hb.ResourceGroup,
    shifted_vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    shift_back_chain: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Lifts over shifted vcf of control region and combines it with the rest of the chrM calls.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('genotype_mito', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        picard LiftoverVcf \
        I={shifted_vcf.vcf} \
        O={j.lifted_vcf} \
        R={reference.fasta} \
        CHAIN={shift_back_chain} \
        REJECT={j.rejected_vcf}

        picard MergeVcfs \
        I=~{vcf.vcf} \
        I=~{j.lifted_vcf} \
        O=~{j.rejected_vcf}.merged.vcf
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_cram, str(output_cram_path.with_suffix('')))

    return j


def genotype_mito(
    b,
    cram_path: CramPath,
    shifted_cram_path: CramPath,
    output_cram_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Genotype mito genome using mutect2

    """
    # If we can resuse existing output, return the existing cram
    if can_reuse(
        [
            output_cram_path,
        ],
        overwrite,
    ):
        return None, b.read_input_group(
            **{
                'cram': str(output_cram_path),
                'cram.crai': str(output_cram_path) + '.crai',
            }
        )

    job_attrs = job_attrs or {}
    j = b.new_job('genotype_mito', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    # reference = fasta_res_group(b)
    # j.declare_resource_group(
    #     output_cram={
    #         'cram': '{root}.cram',
    #         'cram.crai': '{root}.cram.crai',
    #     }
    # )

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
            -O {j.output_cram.cram}

        ls -l
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_cram, str(output_cram_path.with_suffix('')))

    return j, j.output_cram
