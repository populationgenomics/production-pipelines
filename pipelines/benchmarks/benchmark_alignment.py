#!/usr/bin/env python3

"""
Benchmarking different aligment setups.
"""

from enum import Enum
import click
import logging

from cpg_utils.config import update_dict, get_config
from cpg_utils.hail_batch import image_path, fasta_res_group

from cpg_pipes import benchmark, Namespace
from cpg_pipes.jobs.align import Aligner, MarkDupTool, align
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.types import FastqPair, CramPath, SequencingType, FastqPairs
from cpg_pipes.pipeline import (
    stage,
    SampleStage,
    StageInput,
    StageOutput,
)
from cpg_pipes.targets import Sample

logger = logging.getLogger(__file__)


DATASET = 'fewgenomes'
NAMESPACE = Namespace.MAIN


class InputsType(Enum):
    """
    Scope of benchmark
    """
    TOY = 'toy'
    FULL = 'full'
    FULL_SUBSET = 'full_subset'


INPUTS_TYPE = InputsType.FULL


@stage
class SubsetAlignmentInput(SampleStage):
    """
    Subset FASTQ to target number of reads
    """
    def expected_outputs(self, sample: Sample):
        """
        Expected to generate a pair of FASTQs
        """
        basepath = f'{benchmark.BENCHMARK_BUCKET}/outputs/{sample.id}/subset'
        return {
            'r1': f'{basepath}/R1.fastq.gz',
            'r2': f'{basepath}/R2.fastq.gz',
        }

    NA12878_read_pairs = 787265109
    FQ_LINES_PER_READ = 4
    SUBSET_FRACTION = 4

    def _subset_fastq(self, alignment_input, sample):
        fqs1, fqs2 = alignment_input.as_fq_inputs(self.b)

        # We want to take 1/4 of all reads in NA12878. Total read count is 1574530218
        # (787265109 pairs), so subsetting to 196816277 pairs (196816277*4=787265108
        # lines from each file)
        lines = self.NA12878_read_pairs * self.FQ_LINES_PER_READ / self.SUBSET_FRACTION

        j1 = self.b.new_job('Subset FQ1', dict(sample=sample.run_id))
        j1.storage('100G')
        j1.command(f'gunzip -c {fqs1[0]} | head -n{lines} | gzip -c > {j1.out_fq}')
        self.b.write_output(j1.out_fq, str(self.expected_outputs(sample)['r1']))

        j2 = self.b.new_job('Subset FQ2')
        j2.storage('100G')
        j2.command(f'gunzip -c {fqs2[0]} | head -n{lines} | gzip -c > {j2.out_fq}')
        self.b.write_output(j2.out_fq, str(self.expected_outputs(sample)['r2']))

        return self.make_outputs(
            sample,
            data={'r1': j1.out_fq, 'r2': j2.out_fq},
            jobs=[j1, j2],
        )

    def _subset_cram(self, cram: CramPath, sample: Sample):
        j = self.b.new_job('Subset CRAM')
        j.image(image_path('samtools'))
        j.storage('100G')
        reference = fasta_res_group(self.b)
        cram_group = cram.resource_group(self.b)

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
        self.b.write_output(
            j.output_cram, str(self.expected_outputs(sample)).replace('.cram', '')
        )

        return self.make_outputs(
            sample,
            data={'cram': j.output_cram},
            jobs=[j],
        )

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput | None:
        """Queue jobs"""
        if 'fastq_input' in sample.meta:
            alignment_input = sample.meta['fastq_input']
            return self._subset_fastq(alignment_input, sample)
        else:
            return self.make_outputs(sample, None)


@stage(
    required_stages=SubsetAlignmentInput
    if (INPUTS_TYPE == InputsType.FULL_SUBSET)
    else None
)
class DifferentResources(SampleStage):
    """
    Benchmark different resources
    """
    def expected_outputs(self, sample: 'Sample'):
        """
        Not exptected to produce anything
        """
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput | None:
        """
        Add align jobs for each combination of resources/tools
        """
        basepath = (
            benchmark.BENCHMARK_BUCKET / f'outputs/{INPUTS_TYPE.value}/{sample.id}'
        )

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetAlignmentInput)
            alignment_input = FastqPairs([FastqPair(d['r1'], d['r2'])])
        else:
            alignment_input = sample.meta['fastq_input']

        jobs = []
        for nthreads in [8, 16, 32]:
            for aligner in [Aligner.DRAGMAP, Aligner.BWA]:
                jobs.extend(
                    align(
                        self.b,
                        alignment_input=alignment_input,
                        sample_name=sample.id,
                        output_path=basepath
                        / f'nomarkdup/{aligner.name}_nthreads{nthreads}.bam',
                        job_attrs=sample.get_job_attrs(),
                        aligner=aligner,
                        markdup_tool=MarkDupTool.NO_MARKDUP,
                        extra_label=f'nomarkdup_fromfastq_{aligner.name}nthreads{nthreads}',
                        requested_nthreads=nthreads,
                    )
                )
        return self.make_outputs(sample, jobs=jobs)


@stage(
    required_stages=SubsetAlignmentInput
    if (INPUTS_TYPE == InputsType.FULL_SUBSET)
    else None
)
class DifferentAlignerSetups(SampleStage):
    """
    Try different setups of tools (alingners and deduplicators)
    """
    def expected_outputs(self, sample: 'Sample'):
        """
        Not exptected to produce anything
        """
        return None

    def queue_jobs(self, sample: 'Sample', inputs: StageInput) -> StageOutput | None:
        """
        Add align jobs for each combination of tools
        """
        basepath = benchmark.BENCHMARK_BUCKET / f'outputs/{sample.id}'

        if INPUTS_TYPE == InputsType.FULL_SUBSET:
            d = inputs.as_dict(sample, stage=SubsetAlignmentInput)
            alignment_input = FastqPairs([FastqPair(d['r1'], d['r2'])])
        else:
            if 'fastq_input' in sample.meta:
                alignment_input = FastqPairs(sample.meta['fastq_input'])
            else:
                alignment_input = FastqPairs(sample.meta['cram_input'])

        jobs = []
        for aligner in [
            Aligner.DRAGMAP,
            Aligner.BWA,
        ]:
            for markdup in [MarkDupTool.PICARD]:
                jobs.extend(
                    align(
                        self.b,
                        alignment_input=alignment_input,
                        sample_name=sample.id,
                        job_attrs=sample.get_job_attrs(),
                        output_path=basepath / f'{aligner.name}-{markdup.name}.bam',
                        aligner=aligner,
                        markdup_tool=markdup,
                        extra_label=f', dedup with {markdup.name}',
                        realignment_shards_num=10
                        if isinstance(alignment_input, CramPath)
                        else 0,
                    )
                )
        return self.make_outputs(sample, jobs=jobs)


@click.command()
def main():
    """
    Entry point.
    """
    update_dict(get_config()['workflow'], {
        'name': 'Benchmark alignment',
        'dataset': DATASET,
        'access_level': 'full',
    })
    pipeline = Pipeline()

    p = pipeline.create_dataset('fewgenomes')
    s = p.add_sample(
        id='PERTHNEURO_CRAM',
        external_id='PERTHNEURO_CRAM',
    )
    s.alignment_input_by_seq_type[SequencingType.GENOME] = benchmark.perth_neuro_cram

    pipeline.set_stages(
        [
            DifferentAlignerSetups,
        ]
    )

    pipeline.run()


if __name__ == '__main__':
    main()
