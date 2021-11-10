#!/usr/bin/env python3
import click
import logging

from cpg_production_pipelines.hailbatch import AlignmentInput
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.align import Aligner, MarkDupTool, align

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
        title='Align CPG12062_19W001482_A0131064_proband',
        smdb_check_existence=False,
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
        markdup_tool=MarkDupTool.BIOBAMBAM,
    )
    # produce_gvcf(
    #     dragen_mode=True,
    # )
    pipe.run()


if __name__ == '__main__':
    main()
