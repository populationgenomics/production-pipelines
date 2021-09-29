#!/usr/bin/env python3

import hailtop.batch as hb
from os.path import join, splitext
import os
from google.cloud import storage
from typing import Optional, List
from dataclasses import dataclass
import logging

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
BWAMEM2_IMAGE = f'{AR_REPO}/alignment:v4'
DRAGMAP_IMAGE = f'{AR_REPO}/dragmap:1.2.1'
PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'

REF_BUCKET = 'gs://cpg-reference/hg38/v1'
REF_FASTA = join(REF_BUCKET, 'Homo_sapiens_assembly38.fasta')

DRAGMAP_INDEX_BUCKET = 'gs://cpg-reference/dragmap/v0'
DRAGMAP_INDEX_FILES = [
    'hash_table.cfg',
    'hash_table.cfg.bin',
    'hash_table.cmp',
    'hash_table_stats.txt',
    'ref_index.bin',
    'reference.bin',
    'repeat_mask.bin',
    'str_table.bin',
]


@dataclass
class AlignmentInput:
    """
    Sort of a union type for possible alignment inputs
    """

    bam_or_cram_path: Optional[str] = None
    index_path: Optional[str] = None
    fqs1: Optional[List[str]] = None
    fqs2: Optional[List[str]] = None


def _test_dragmap(
    b: hb.Batch,
    sample_name: str,
    alignment_input: AlignmentInput,
    enable_sampling: bool = False,
):
    dragmap_index = b.read_input_group(**{
        k.replace('.', '_'): join(DRAGMAP_INDEX_BUCKET, k) for k in DRAGMAP_INDEX_FILES
    })
    
    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )
 
    label = f'{sample_name}-dragmap-{"with-enable-sampling" if enable_sampling else ""}'
    if alignment_input.bam_or_cram_path:
        label += '-from-cram'
    else:
        label += '-from-fq'
    
    output_cram_path = f'gs://cpg-seqr-test-tmp/test-alignment/{label}.cram'
    job_name = f'Test DRAGMAP {label}'
    j = b.new_job(job_name)
    j.image(DRAGMAP_IMAGE)
    j.cpu(32)
    j.storage('300G')

    if alignment_input.bam_or_cram_path:
        cram = b.read_input_group(
            base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
        )
        input_param = f'--bam-input {cram} --ht-reference {reference.base}'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        input_param = f'-1 <(cat {" ".join(files1)}) -2 <(cat {" ".join(files2)})'
    logger.info(f'input_param: {input_param}')

    command = f"""
set -o pipefail
set -ex
(while true; do df -h; pwd; du -sh $(dirname {j.sorted_bam}); sleep 600; done) &

dragen-os -r {dragmap_index} {input_param} \\
--RGID {sample_name} --RGSM {sample_name} \\
{"--enable-sampling=1" if enable_sampling else ""} |
samtools sort -T $(dirname {j.sorted_bam})/samtools-sort-tmp -Obam -o {j.sorted_bam}

df -h; pwd; du -sh $(dirname {j.sorted_bam})
    """
    j.command(command)

    md_j = b.new_job('MarkDuplicates')
    md_j.image(PICARD_IMAGE)
    md_j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'cram.crai': '{root}.cram.crai',
        }
    )
    command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {md_j.output_cram.cram}); sleep 600; done) &

picard MarkDuplicates \\
    I={j.sorted_bam} O=/dev/stdout M={md_j.duplicate_metrics} \\
    TMP_DIR=$(dirname {md_j.output_cram.cram})/picard-tmp \\
    ASSUME_SORT_ORDER=coordinate | \\
samtools view -@30 -T {reference.base} -O cram -o {md_j.output_cram.cram}

samtools index -@2 {md_j.output_cram.cram} {md_j.output_cram["cram.crai"]}

df -h; pwd; du -sh $(dirname {md_j.output_cram.cram})
    """
    md_j.command(command)
    md_j.cpu(2)
    md_j.memory('standard')
    md_j.storage('150G')
    b.write_output(md_j.output_cram, splitext(output_cram_path)[0])
    return md_j, output_cram_path


def _test_bwamem2(
    b: hb.Batch,
    sample_name: str,
    alignment_input: AlignmentInput,
    use_picard: bool = False,
    use_dragmap: bool = False,
):
    bwa_reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
        sa=REF_FASTA + '.sa',
        amb=REF_FASTA + '.amb',
        bwt=REF_FASTA + '.bwt',
        ann=REF_FASTA + '.ann',
        pac=REF_FASTA + '.pac',
        o123=REF_FASTA + '.0123',
        bwa2bit64=REF_FASTA + '.bwt.2bit.64',
    )
    
    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )

    label = f'{sample_name}-bwamem2-{"picard" if use_picard else "biobamba"}' 
    output_cram_path = f'gs://cpg-seqr-test-tmp/test-alignment/{label}.cram'
    job_name = f'Test BWA {label}'
    j = b.new_job(job_name)
    j.image(DRAGMAP_IMAGE if use_dragmap else BWAMEM2_IMAGE)
    total_cpu = 32
    bwa_cpu = total_cpu

    if alignment_input.bam_or_cram_path:
        use_bazam = True
        bazam_cpu = 10
        assert alignment_input.index_path
        assert not alignment_input.fqs1 and not alignment_input.fqs2
        cram = b.read_input_group(
            base=alignment_input.bam_or_cram_path, index=alignment_input.index_path
        )
        r1_param = (
            f'<(bazam -Xmx16g -Dsamjdk.reference_fasta={bwa_reference.base}'
            f' -n{bazam_cpu} -bam {cram.base})'
        )
        r2_param = '-'
    else:
        assert alignment_input.fqs1 and alignment_input.fqs2
        use_bazam = False
        files1 = [b.read_input(f1) for f1 in alignment_input.fqs1]
        files2 = [b.read_input(f1) for f1 in alignment_input.fqs2]
        r1_param = f'<(cat {" ".join(files1)})'
        r2_param = f'<(cat {" ".join(files2)})'
        logger.info(f'r1_param: {r1_param}')
        logger.info(f'r2_param: {r2_param}')

    j.cpu(total_cpu)
    j.memory('standard')
    j.storage('300G')

    rg_line = f'@RG\\tID:{sample_name}\\tSM:{sample_name}'
    # BWA command options:
    # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
    # -p     smart pairing (ignoring in2.fq)
    # -t16   threads
    # -Y     use soft clipping for supplementary alignments
    # -R     read group header line such as '@RG\tID:foo\tSM:bar'
    if not use_picard:
        j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'crai': '{root}.crai',
            }
        )
        command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.output_cram.cram}); sleep 600; done) &

bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
    -R '{rg_line}' {reference.base} {r1_param} {r2_param} | \\
bamsormadup inputformat=sam threads=10 SO=coordinate \\
    M={j.duplicate_metrics} outputformat=sam \\
    tmpfile=$(dirname {j.output_cram.cram})/bamsormadup-tmp | \\
samtools view -T {reference.base} -O cram -o {j.output_cram.cram}

samtools index -@{total_cpu} {j.output_cram.cram} {j.output_cram.crai}

df -h; pwd; du -sh $(dirname {j.output_cram.cram})
        """
        j.command(command)
        b.write_output(j.output_cram, splitext(output_cram_path)[0])
        return j, output_cram_path

    else:
        command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {j.sorted_bam}); sleep 600; done) &

bwa-mem2 mem -K 100000000 {'-p' if use_bazam else ''} -t{bwa_cpu} -Y \\
        -R '{rg_line}' {reference.base} {r1_param} {r2_param} > {j.sorted_bam}

df -h; pwd; du -sh $(dirname {j.sorted_bam})
        """
        j.command(command)

        md_j = b.new_job('MarkDuplicates')
        md_j.image(PICARD_IMAGE)
        md_j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'crai': '{root}.crai',
            }
        )
        command = f"""
set -o pipefail
set -ex

(while true; do df -h; pwd; du -sh $(dirname {md_j.output_cram.cram}); sleep 600; done) &

picard MarkDuplicates \\
        I={j.sorted_bam} O=/dev/stdout M={md_j.duplicate_metrics} \\
        TMP_DIR=$(dirname {md_j.output_cram.cram})/picard-tmp \\
        ASSUME_SORT_ORDER=coordinate | \\
    samtools view -@30 -T {reference.base} -O cram -o {md_j.output_cram.cram}

samtools index -@2 {md_j.output_cram.cram} {md_j.output_cram.crai}

df -h; pwd; du -sh $(dirname {md_j.output_cram.cram})
        """
        md_j.command(command)
        md_j.cpu(2)
        md_j.memory('standard')
        md_j.storage('150G')
        b.write_output(md_j.output_cram, splitext(output_cram_path)[0])
        return md_j, output_cram_path


def file_exists(path: str) -> bool:
    """
    Check if the object exists, where the object can be:
        * local file
        * local directory
        * Google Storage object
        * Google Storage URL representing a *.mt or *.ht Hail data,
          in which case it will check for the existence of a
          *.mt/_SUCCESS or *.ht/_SUCCESS file.
    :param path: path to the file/directory/object/mt/ht
    :return: True if the object exists
    """
    if path.startswith('gs://'):
        bucket = path.replace('gs://', '').split('/')[0]
        path = path.replace('gs://', '').split('/', maxsplit=1)[1]
        path = path.rstrip('/')  # '.mt/' -> '.mt'
        if any(path.endswith(f'.{suf}') for suf in ['mt', 'ht']):
            path = os.path.join(path, '_SUCCESS')
        gs = storage.Client()
        return gs.get_bucket(bucket).get_blob(path)
    return os.path.exists(path)


billing_project = os.getenv('HAIL_BILLING_PROJECT') or 'seqr'
hail_bucket = os.environ.get('HAIL_BUCKET', 'cpg-seqr-test-tmp')
print(
    f'Starting hail Batch with the project {billing_project}, ' f'bucket {hail_bucket}'
)
backend = hb.ServiceBackend(
    billing_project=billing_project,
    bucket=hail_bucket.replace('gs://', ''),
)
b = hb.Batch(backend=backend, name='Benchmark alignment')


# --enable-sampling arg (=1) Automatically detect
#                            paired-end parameters by
#                            running a sample through the
#                            mapper/aligner
for (use_dragmap, enable_sampling, bam_input) in [
    (True, False, False),
    (True, True, False),
    # (True, True, True),
]:
    if not bam_input:
        alignment_input = AlignmentInput(
            fqs1=['gs://cpg-perth-neuro-main-upload/seqr_transfers/HYH7YCCXY_5_190317_FD03043182_Homo-sapiens__R_190312_GINRAV_DNA_M001_R1.fastq.gz'],
            fqs2=['gs://cpg-perth-neuro-main-upload/seqr_transfers/HYH7YCCXY_5_190317_FD03043182_Homo-sapiens__R_190312_GINRAV_DNA_M001_R2.fastq.gz'],
        )
    else:
        alignment_input = None
    if alignment_input:
        logger.info(f'Submitting {alignment_input}')
        _test_dragmap(
            b=b,
            sample_name='HYH7YCCXY',
            alignment_input=alignment_input,
            enable_sampling=enable_sampling,
        )

b.run(wait=False)
