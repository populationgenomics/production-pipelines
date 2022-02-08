#!/usr/bin/env python3
import click
import logging

from cpg_pipes import images
from cpg_pipes.hb.inputs import AlignmentInput
from cpg_pipes.pipeline import Pipeline
from cpg_pipes.jobs.align import Aligner, MarkDupTool, align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


PROJECT = 'seqr'
NAMESPACE = 'main'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='seqr',
        name='seqr_align_CPG12062_19W001482_A0131064_proband',
        output_version='v0',
        namespace=NAMESPACE,
        title='CPG12062_19W001482_A0131064_proband',
        check_smdb_seq_existence=False,
        keep_scratch=True,
    )

    dir_path = 'gs://cpg-acute-care-main-upload/cpg_acute_positives_20211003_213917'
    fq_input = AlignmentInput(
        fqs1=[
            f'{dir_path}/191129_A00692_0037_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L001_R1.fastq.gz',
            f'{dir_path}/191201_A00692_0038_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L001_R1.fastq.gz',
            f'{dir_path}/191129_A00692_0037_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            f'{dir_path}/191201_A00692_0038_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L002_R1.fastq.gz',
        ],
        fqs2=[
            f'{dir_path}/191129_A00692_0037_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L001_R2.fastq.gz',
            f'{dir_path}/191201_A00692_0038_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L001_R2.fastq.gz',
            f'{dir_path}/191129_A00692_0037_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L002_R2.fastq.gz',
            f'{dir_path}/191201_A00692_0038_TL1911296_19W001482-FAM000327_MAN-20191129_NEXTERAFLEXWGS_L002_R2.fastq.gz',
        ]
    )
    align(
        pipe.b,
        alignment_input=fq_input,
        sample_name='CPG12062_19W001482_A0131064_proband',
        output_path=f'gs//cpg-seqr-test/test/CPG12062_19W001482_A0131064.bam',
        project_name=PROJECT,
        aligner=Aligner.BWA,
        markdup_tool=MarkDupTool.PICARD,
    )
    
    merged_bam = pipe.b.read_input('gs://cpg-seqr-main-tmp/seqr_align_CPG12062_19W001482_A0131064_proband/v0/hail/batch/3b4ddd/1/sorted_bam')
    # picard.markdup(
    #     pipe.b,
    #     merged_bam),
    #     sample_name='CPG12062_19W001482_A0131064_proband',
    #     project_name=PROJECT,
    #     overwrite=False,
    # )
    j = pipe.b.new_job('Index bam')
    j.declare_resource_group(
        output={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bam.bai',
        }
    )
    j.storage('300G')
    j.cpu(8)
    j.image(images.SAMTOOLS_PICARD_IMAGE)
    j.command(f"""\
    mv {merged_bam} {j.output.bam}
    samtools index -@7 {j.output.bam} {j.output['bam.bai']}
    """)
    pipe.b.write_output(j.output, 'gs://cpg-seqr-main-tmp/seqr_align_CPG12062_19W001482_A0131064_proband-sorted')

    pipe.submit_batch()


if __name__ == '__main__':
    main()
