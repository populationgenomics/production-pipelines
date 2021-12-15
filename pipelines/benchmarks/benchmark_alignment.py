#!/usr/bin/env python3
from enum import Enum

import click
import logging

from cpg_pipes import benchmark
from cpg_pipes.hailbatch import AlignmentInput
from cpg_pipes.pipeline import Pipeline, stage, SampleStage, StageInput, \
    StageOutput, Sample
from cpg_pipes.jobs.align import Aligner, MarkDupTool, align

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


PROJECT = 'fewgenomes'
NAMESPACE = 'main'


class InputsType(Enum):
    TOY = 1
    FULL = 2
    FULL_SUBSET = 3


INPUTS_TYPE = InputsType.TOY


@stage
class SubsetFastq(SampleStage):
    def expected_result(self, sample: 'Sample'):
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}/subset'
        return {
            'r1': f'{basepath}/R1.fastq.gz',
            'r2': f'{basepath}/R2.fastq.gz',
        }

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:

        alignment_input = sample.meta['fastq_input']

        fqs1, fqs2 = alignment_input.as_fq_inputs(self.pipe.b)

        # We want to take 1/4 of all reads. Total read count is 1574530218 
        # (787265109 pairs), so subsetting to 196816277 pairs (196816277*4=787265108 
        # lines from each file)
        lines = 787265108

        j1 = self.pipe.b.new_job('Subset FQ1', dict(sample=sample.id))
        j1.storage('100G')
        j1.command(f'gunzip -c {fqs1[0]} | head -n{lines} | gzip -c > {j1.out_fq}')
        self.pipe.b.write_output(j1.out_fq, self.expected_result(sample)['r1'])

        j2 = self.pipe.b.new_job('Subset FQ2')
        j2.storage('100G')
        j2.command(f'gunzip -c {fqs2[0]} | head -n{lines} | gzip -c > {j2.out_fq}')
        self.pipe.b.write_output(j2.out_fq, self.expected_result(sample)['r2'])
        
        return self.make_outputs(
            sample, 
            data={'r1': j1.out_fq, 'r2': j2.out_fq},
            jobs=[j1, j2],
        )


@stage(requires_stages=SubsetFastq if INPUTS_TYPE == InputsType.FULL_SUBSET else None)
class DifferentResources(SampleStage):
    def expected_result(self, sample: 'Sample'):
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}'

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetFastq)
            alignment_input = AlignmentInput(fqs1=[d['r1']], fqs2=[d['r2']])
        else:
            alignment_input = sample.meta['fastq_input']

        jobs = []
        for nthreads in [8, 16, 32]:
            for aligner in [Aligner.DRAGMAP, Aligner.BWA]:
                jobs.append(align(
                    self.pipe.b,
                    alignment_input=alignment_input,
                    sample_name=sample.id,
                    output_path=f'{basepath}/nomarkdup/{aligner.name}_nthreads{nthreads}.bam',
                    project_name=PROJECT,
                    aligner=aligner,
                    markdup_tool=MarkDupTool.NO_MARKDUP,
                    extra_label=f'nomarkdup_fromfastq_{aligner.name}nthreads{nthreads}',
                    depends_on=inputs.get_jobs(),
                    nthreads=nthreads,
                ))
        return self.make_outputs(sample, jobs=jobs)


@stage(requires_stages=SubsetFastq if INPUTS_TYPE == InputsType.FULL_SUBSET else None)
class DifferentAlignerSetups(SampleStage):
    def expected_result(self, sample: 'Sample'):
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}'

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetFastq)
            alignment_input = AlignmentInput(fqs1=[d['r1']], fqs2=[d['r2']])
        else:
            alignment_input = sample.meta['fastq_input']

        jobs = []
        for aligner in [Aligner.DRAGMAP, Aligner.BWA, Aligner.BWAMEM2]:
            for markdup in [MarkDupTool.PICARD, MarkDupTool.BIOBAMBAM]:
                jobs.append(align(
                    self.pipe.b,
                    alignment_input=alignment_input,
                    sample_name=sample.id,
                    project_name=PROJECT,
                    output_path=f'{basepath}/{aligner.name}-{markdup.name}.bam',
                    aligner=aligner,
                    markdup_tool=markdup,
                    extra_label=f'{aligner.name}-{markdup.name}',
                    depends_on=inputs.get_jobs(),
                ))
        return self.make_outputs(sample, jobs=jobs)


@click.command()
def main():
    pipeline = Pipeline(
        name='benchmark_alignment',
        title='Benchmark alignment',
        analysis_project=PROJECT,
        output_version='v0',
        namespace=NAMESPACE,
        check_smdb_seq_existence=False,
        keep_scratch=True,
    )

    p = pipeline.add_project('fewgenomes')
    p.add_sample(
        id='NA12340',
        external_id='NA12340',
        fastq_input=benchmark.na12878fq
    )

    pipeline.set_stages([
        SubsetFastq,
        DifferentAlignerSetups,
    ])

    pipeline.submit_batch()


if __name__ == '__main__':
    main()
