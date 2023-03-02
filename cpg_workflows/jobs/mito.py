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


def mito_realign(
    b,
    sample: Sample,
    cram_path: CramPath,
    mt_aligned_cram: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Re-align reads to mito genome
    """
    if can_reuse(
        [
            cram_path,
        ],
        overwrite,
    ):
        return None

    job_attrs = job_attrs or {}
    j = b.new_job('mito_realign', job_attrs)
    j.image(image_path('bwamem2'))

    reference = fasta_res_group(b)

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

    shifted_mt_ref = b.read_input_group(
        dict='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.shifted_by_8000_bases.chrM.dict',
        fasta='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.shifted_by_8000_bases.chrM.fasta',
        amb='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb',
        ann='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann',
        bwt='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt',
        fai='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai',
        pac='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac',
        sa='gs://cpg-common-main/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa',
    )

    # As we are only accessing a tiny chunk of the cram
    # accessing the cram via cloudfuse is faster than localising the full cram
    bucket = cram_path.path.drive
    print(f'bucket = {bucket}')
    bucket_mount_path = to_path('/bucket')
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)
    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])
    assert cram_path.index_path  # keep mypy happy as index_path is optional
    mounted_cram_index_path = bucket_mount_path / '/'.join(
        cram_path.index_path.parts[2:]
    )

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    nthreads = res.get_nthreads()

    cmd = dedent(
        f"""\
        bazam -Xmx16g -Dsamjdk.reference_fasta={reference.base} \
            -n{min(nthreads, 6)} {mounted_cram_path} -L chrM | \
         bwa mem -p {mt_ref.fasta}  - | \
         samtools view -bSu - | \
         samtools sort -o out.bam
        """
    )

    j.command(command(cmd))
    b.write_output(j.mt_aligned_cram, str(mt_aligned_cram))

    return j
