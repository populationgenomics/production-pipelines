#!/usr/bin/env python3

"""
Benchmarking different aligment setups.
"""

from enum import Enum
import click
import logging

from cpg_pipes import benchmark, images
from cpg_pipes.jobs.align import Aligner, MarkDupTool, align
from cpg_pipes.pipeline.analysis import FastqPair, CramPath
from cpg_pipes.pipeline.pipeline import stage, Pipeline
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput
from cpg_pipes.ref_data import REF_D

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


DATASET = 'fewgenomes'
NAMESPACE = 'main'


class InputsType(Enum):
    TOY = 'toy'
    FULL = 'full'
    FULL_SUBSET = 'full_subset'


INPUTS_TYPE = InputsType.FULL


@stage
class SubsetAlignmentInput(SampleStage):
    def expected_result(self, sample: Sample):
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}/subset'
        return {
            'r1': f'{basepath}/R1.fastq.gz',
            'r2': f'{basepath}/R2.fastq.gz',
        }
    
    NA12878_read_pairs = 787265109
    FQ_LINES_PER_READ = 4
    SUBSET_FRACTION = 4
    
    def _subset_fastq(self, alignment_input, sample):
        fqs1, fqs2 = alignment_input.as_fq_inputs(self.pipe.b)

        # We want to take 1/4 of all reads in NA12878. Total read count is 1574530218 
        # (787265109 pairs), so subsetting to 196816277 pairs (196816277*4=787265108 
        # lines from each file)
        lines = self.NA12878_read_pairs * self.FQ_LINES_PER_READ / self.SUBSET_FRACTION

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
    
    def _subset_cram(self, cram: CramPath, sample: Sample):
        j = self.pipe.b.new_job('Subset CRAM')
        j.image(images.BIOINFO_IMAGE)
        j.storage('100G')
        reference = self.pipe.b.read_input_group(**REF_D)
        cram_group = cram.resource_group(self.pipe.b)

        j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            }
        )
        j.command(
            f'samtools view {cram_group.cram} -@31 --subsample {1/self.SUBSET_FRACTION} '
            f'-T {reference.base} -Ocram -o {j.output_cram.cram_path}\n'
            f'samtools index -@31 {j.output_cram.cram_path} {j.output_cram["cram.crai"]}'
        )
        self.pipe.b.write_output(
            j.output_cram, 
            self.expected_result(sample).replace('.cram', '')
        )
        
        return self.make_outputs(
            sample, 
            data={'cram': j.output_cram},
            jobs=[j],
        )
    
    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        if 'fastq_input' in sample.meta:
            alignment_input = sample.meta['fastq_input']
            return self._subset_fastq(alignment_input, sample)
        else:
            return self.make_outputs(sample, None)


@stage(required_stages=SubsetAlignmentInput if (INPUTS_TYPE == InputsType.FULL_SUBSET) else None)
class DifferentResources(SampleStage):
    def expected_result(self, sample: 'Sample'):
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{INPUTS_TYPE.value}/{sample.id}'

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetAlignmentInput)
            alignment_input = [FastqPair(d['r1'], d['r2'])]
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
                    dataset_name=DATASET,
                    aligner=aligner,
                    markdup_tool=MarkDupTool.NO_MARKDUP,
                    extra_label=f'nomarkdup_fromfastq_{aligner.name}nthreads{nthreads}',
                    depends_on=inputs.get_jobs(),
                    requested_nthreads=nthreads,
                ))
        return self.make_outputs(sample, jobs=jobs)


@stage(required_stages=SubsetAlignmentInput if (INPUTS_TYPE == InputsType.FULL_SUBSET) else None)
class DifferentAlignerSetups(SampleStage):
    def expected_result(self, sample: 'Sample'):
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput:
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}'

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetAlignmentInput)
            alignment_input = [FastqPair(d['r1'], ['r2'])]
        else:
            if 'fastq_input' in sample.meta:
                alignment_input = sample.meta['fastq_input']
            else:
                alignment_input = sample.meta['cram_input']

        jobs = []
        for aligner in [
            Aligner.DRAGMAP, 
            Aligner.BWA,
            # Aligner.BWAMEM2
        ]:
            for markdup in [MarkDupTool.PICARD]:
                jobs.append(align(
                    self.pipe.b,
                    alignment_input=alignment_input,
                    sample_name=sample.id,
                    dataset_name=DATASET,
                    output_path=f'{basepath}/{aligner.name}-{markdup.name}.bam',
                    aligner=aligner,
                    markdup_tool=markdup,
                    extra_label=f', dedup with {markdup.name}',
                    depends_on=inputs.get_jobs(),
                    number_of_shards_for_realignment=10 if 
                    isinstance(alignment_input, CramPath) else 0,
                ))
        return self.make_outputs(sample, jobs=jobs)


@click.command()
def main():
    pipeline = Pipeline(
        name='benchmark_alignment',
        description='Benchmark alignment',
        analysis_dataset=DATASET,
        output_version='v0',
        namespace=NAMESPACE,
        check_smdb_seq=False,
        keep_scratch=True,
    )

    p = pipeline.add_dataset('fewgenomes')
    p.add_sample(
        id='PERTHNEURO_CRAM',
        external_id='PERTHNEURO_CRAM',
        cram_input=benchmark.perth_neuro_cram,
    )

    pipeline.set_stages([
        DifferentAlignerSetups,
    ])

    pipeline.submit_batch()


if __name__ == '__main__':
    main()
