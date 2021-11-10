#!/usr/bin/env python3
import click
import logging

from cpg_production_pipelines import benchmark
from cpg_production_pipelines.hailbatch import AlignmentInput
from cpg_production_pipelines.pipeline import Pipeline
from cpg_production_pipelines.jobs.align import Aligner, MarkDupTool, align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


PROJECT = 'fewgenomes'
NAMESPACE = 'main'


def _subset_fq(pipe, sn: str, inp: AlignmentInput):
    # We want to take 1/4 of all reads. Total read count is 1574530218 
    # (787265109 pairs), so subsetting to 196816277 pairs (196816277*4=787265108 
    # lines from each file)
    lines = 787265108
    fqs1, fqs2 = inp.get_fq_inputs(pipe.b)
    j1 = pipe.b.new_job('Subset FQ1')
    j1.storage('100G')
    j1.command(f'gunzip -c {fqs1[0]} | head -n{lines} | gzip -c > {j1.out_fq}')
    pipe.b.write_output(j1.out_fq, f'{benchmark.BENCHMARK_BUCKET}/{sn}/r1.fq.gz')
    j2 = pipe.b.new_job('Subset FQ2')
    j2.storage('100G')
    j2.command(f'gunzip -c {fqs2[0]} | head -n{lines} | gzip -c > {j2.out_fq}')
    pipe.b.write_output(j2.out_fq, f'{benchmark.BENCHMARK_BUCKET}/{sn}/r2.fq.gz')
    inp = AlignmentInput(
        fqs1=[f'{benchmark.BENCHMARK_BUCKET}/{sn}/r1.fq.gz'], 
        fqs2=[f'{benchmark.BENCHMARK_BUCKET}/{sn}/r2.fq.gz'],
    )
    deps = [j1, j2]
    return inp, deps


@click.command()
def main():
    pipe = Pipeline(
        analysis_project='fewgenomes',
        name='benchmark_alignment',
        output_version='v0',
        namespace=NAMESPACE,
        title='Benchmark alignment',
        smdb_check_existence=False,
        keep_scratch=True,
    )

    fq_inputs = {
        'NA12878': benchmark.na12878fq
    }
    
    for sample_name, inp in fq_inputs.items():
        # deps, inp = _subset(pipe, sn, inp)
        deps = []

        for ncpu in [8, 16, 32]:
            align(
                pipe.b,
                alignment_input=inp,
                sample_name=sample_name,
                output_path=f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}/nomarkdup/dragmap.bam',
                project_name=PROJECT,
                aligner=Aligner.DRAGMAP,
                markdup_tool=MarkDupTool.NO_MARKDUP,
                extra_label=f'nomarkdup_fromfastq_ncpu{ncpu}',
                depends_on=deps,
                ncpu=ncpu,
                highmem=True if ncpu < 32 else False,
                storage_gb=200,
            )
            align(
                pipe.b,
                alignment_input=inp,
                sample_name=sample_name,
                output_path=f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample_name}/nomarkdup/bwamem.bam',
                project_name=PROJECT,
                aligner=Aligner.BWA,
                markdup_tool=MarkDupTool.NO_MARKDUP,
                extra_label=f'nomarkdup_fromfastq_ncpu{ncpu}',
                depends_on=deps,
                ncpu=ncpu,
                highmem=True if ncpu < 32 else False,
                storage_gb=200,
            )

        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     sample_name=sn,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-picard.bam',
        #     project_name='Benchmark',
        #     aligner=Aligner.DRAGMAP,
        #     markdup_tool=MarkDupTool.PICARD,
        #     extra_label='picard',
        # )
        # 
        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-picard.bam',
        #     sample_name=sn,
        #     project_name='Benchmark',
        #     aligner=Aligner.BWA,
        #     markdup_tool=MarkDupTool.PICARD,
        #     extra_label='picard',
        # )
        # 
        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-picard.bam',
        #     sample_name=sn,
        #     project_name='Benchmark',
        #     aligner=Aligner.BWAMEM2,
        #     markdup_tool=MarkDupTool.PICARD,
        #     extra_label='picard',
        # )
        # 
        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/dragmap-biobambam.bam',
        #     sample_name=sn,
        #     project_name='Benchmark',
        #     aligner=Aligner.DRAGMAP,
        #     markdup_tool=MarkDupTool.BIOBAMBAM,
        #     extra_label='biobambam',
        # )
        # 
        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/bwa-biobambam.bam',
        #     sample_name=sn,
        #     project_name='Benchmark',
        #     aligner=Aligner.BWA,
        #     markdup_tool=MarkDupTool.BIOBAMBAM,
        #     extra_label='biobambam',
        # )
        # 
        # align(
        #     pipe.b,
        #     alignment_input=inp,
        #     output_path=f'{BENCHMARK_BUCKET}/{sn}/bwamem2-biobambam.bam',
        #     sample_name=sn,
        #     project_name='Benchmark',
        #     aligner=Aligner.BWAMEM2,
        #     markdup_tool=MarkDupTool.BIOBAMBAM,
        #     extra_label='biobambam',
        # )

        # produce_gvcf(
        #     dragen_mode=True,
        # )
    
    pipe.run()


if __name__ == '__main__':
    main()
