"""
Create Hail Batch jobs to run STRipy
"""

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
    mito_subset_cram: str,
    mito_subset_crai: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Extract read pairs with at least one read mapping to chrM

    This uses the same extraction command as the broad pipeline
    as I am not clear the specifics *may* be important to how
    NUMT mapping reads are handled
    """

    job_attrs = job_attrs or {}
    j = b.new_job('subset_cram_to_chrM', job_attrs)
    j.image(image_path('gatk'))

    reference = fasta_res_group(b)

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

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
            -O {j.mito_subset_cram}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.mito_subset_cram, str(mito_subset_cram))
    b.write_output(j.mito_subset_crai, str(mito_subset_crai))

    return j


def mito_realign(
    b,
    cram_path: CramPath,
    mito_aligned_cram: str,
    mito_aligned_crai: str,
    shifted: bool,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Re-align reads to mito genome
    """
    # if can_reuse(
    #     [
    #         cram_path,
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
            dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.shifted_by_8000_bases.chrM.dict',
            fasta='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.shifted_by_8000_bases.chrM.fasta',
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
        CRAM=$BATCH_TMPDIR/{cram_path.path.name}
        CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path.path)} $CRAM
        retry_gs_cp {str(cram_path.index_path)} $CRAI

        bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
            -n{min(nthreads, 6)} -bam $CRAM -L chrM | \
         bwa mem -K 100000000 -p -v 3 -t 2 -Y {mt_ref.fasta} - | \
         samtools view -bSu - | \
        samtools sort -o {j.raw_cram}
        """
    )

    j.command(command(cmd, define_retry_function=True))

    mkdup_j = picard.markdup(
        b=b,
        sorted_bam=j.raw_cram,
        output_path=mito_aligned_cram,
        out_markdup_metrics_path=mito_aligned_cram + '.markduplicates-metrics')

    return [j, mkdup_j]
