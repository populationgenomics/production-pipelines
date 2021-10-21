#!/usr/bin/env python3
import click
import logging
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.fastqc import fastqc
from cpg_production_pipelines.utils import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


INPUT_PROJECT = 'fewgenomes'
NAMESPACE = 'main'
BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test/benchmark'
INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs'
TOY_INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs/test'
RESULTS_BUCKET = f'{BENCHMARK_BUCKET}/fastqc'


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='run_qc',
        output_version='v0',
        namespace=NAMESPACE,
        title='Run QC',
        smdb_check_existence=False,
    )
    # This samples crashes with BWA-MEM2 after 30min
    cpg13326 = AlignmentInput(
        fqs1=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r1.fastq.gz'],
        fqs2=[f'{BENCHMARK_BUCKET}/inputs/D18-1351_r2.fastq.gz'],
    )
    
    perth_neuro2 = AlignmentInput(
        fqs1=[f'{INPUTS_BUCKET}/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R1.fastq.gz'],
        fqs2=[f'{INPUTS_BUCKET}/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R2.fastq.gz'],
    )
    
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

    # for sample_name, inp in giab_inputs.items():
    # for sample_name, inp in tiny_inputs.items():
    for sample_name, inp in {'PERTH_NEURO2': perth_neuro2}.items():
        fastqc(
            pipe.b,
            results_bucket=RESULTS_BUCKET,
            alignment_input=inp,
            sample_name=sample_name,
            project_name='FastQC',
        )
    
    pipe.run()


if __name__ == '__main__':
    main()
