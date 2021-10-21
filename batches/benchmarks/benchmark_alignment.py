#!/usr/bin/env python3
import click
import logging
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.align import AlignmentInput, Aligner, MarkDupTool, \
    align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_PROJECT = 'fewgenomes'
NAMESPACE = 'main'
BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test/benchmark'
TOY_INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs/test'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='benchmark_alignment',
        output_version='v0',
        namespace=NAMESPACE,
        title='Benchmark alignment',
        smdb_check_existence=False,
    )
    # This samples crashes with BWA-MEM2 after 30min
    # cpg13326 = AlignmentInput(
    #     fqs1=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r1.fastq.gz'],
    #     fqs2=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r2.fastq.gz'],
    # )
    # This set is 50MB each
    tiny_fq = AlignmentInput(
        fqs1=[f'{TOY_INPUTS_BUCKET}/2-699835.L001.R1.n40000.fastq.gz'],
        fqs2=[f'{TOY_INPUTS_BUCKET}/2-699835.L001.R2.n40000.fastq.gz'],
    )
    tiny_cram_path = f'{TOY_INPUTS_BUCKET}/NA12878-chr21-tiny.cram'
    tiny_cram = AlignmentInput(
        bam_or_cram_path=tiny_cram_path,
        index_path=tiny_cram_path + '.crai',
    )
    tiny_inputs = {
        'TINY_FQ': tiny_fq,
        'TINY_CRAM': tiny_cram,
    }
    
    giab_path = 'gs://cpg-reference/validation/giab/cram'
    giab_inputs = {
        sn: AlignmentInput(
            bam_or_cram_path=f'{giab_path}/{sn}.cram',
            index_path=f'{giab_path}/{sn}.cram.crai',
        )
        for sn in ['NA12878', 'NA12891', 'NA12892']
    }

    # for sn, inp in giab_inputs.items():
    for sn, inp in tiny_inputs.items():
        align(
            pipe.b,
            alignment_input=inp,
            sample_name=sn,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-picard.bam',
            project_name='Benchmark',
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.PICARD,
            extra_label='picard',
        )

        align(
            pipe.b,
            alignment_input=inp,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-picard.bam',
            sample_name=sn,
            project_name='Benchmark',
            aligner=Aligner.BWA,
            markdup_tool=MarkDupTool.PICARD,
            extra_label='picard',
        )
    
        align(
            pipe.b,
            alignment_input=inp,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-picard.bam',
            sample_name=sn,
            project_name='Benchmark',
            aligner=Aligner.BWAMEM2,
            markdup_tool=MarkDupTool.PICARD,
            extra_label='picard',
        )
    
        align(
            pipe.b,
            alignment_input=inp,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-biobambam.bam',
            sample_name=sn,
            project_name='Benchmark',
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.BIOBAMBAM,
            extra_label='biobambam',
        )
    
        align(
            pipe.b,
            alignment_input=inp,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-biobambam.bam',
            sample_name=sn,
            project_name='Benchmark',
            aligner=Aligner.BWA,
            markdup_tool=MarkDupTool.BIOBAMBAM,
            extra_label='biobambam',
        )
    
        align(
            pipe.b,
            alignment_input=inp,
            output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-biobambam.bam',
            sample_name=sn,
            project_name='Benchmark',
            aligner=Aligner.BWAMEM2,
            markdup_tool=MarkDupTool.BIOBAMBAM,
            extra_label='biobambam',
        )

        # produce_gvcf(
        #     dragen_mode=True,
        # )
    
    pipe.run()


if __name__ == '__main__':
    main()
