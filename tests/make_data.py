#!/usr/bin/env python3

"""
Generate data for unit tests
"""
from typing import cast

from cpg_pipes import images, Namespace
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs.align import extract_fastq
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.types import SequencingType, CramPath

try:
    from .utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID, setup_env, SUBSET_FQ_BY_SID, FULL_GVCF_BY_SID, FULL_CRAM_BY_SID, SUBSET_CRAM_BY_SID
except ImportError:
    from utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID, setup_env, SUBSET_FQ_BY_SID, FULL_GVCF_BY_SID, FULL_CRAM_BY_SID, SUBSET_CRAM_BY_SID  # type: ignore


def main():
    """
    Generate data for unit tests
    """
    setup_env()
    pipeline = create_pipeline(
        name='make_test_data',
        description='Make test data',
        analysis_dataset=DATASET,
        namespace=Namespace.TEST,
    )
    make_subset_crams(pipeline)
    make_gvcfs_for_joint_calling(pipeline)


def make_subset_crams(pipeline: Pipeline):
    """
    Make toy CRAMs that span entire genome, not just chr20.
    1. Take WGS CRAMs,
    2. Randomly create regions.
    3. Reduce the coverage to ~5x,
    4. Convert to fastq pairs.
    """
    b = pipeline.b
    d = pipeline.cohort.create_dataset(DATASET)
    samples = [d.add_sample(
        sid, external_id=sid, alignment_input=CramPath(FULL_CRAM_BY_SID[sid])
    ) for sid in SAMPLES]

    refs = pipeline.refs

    intervals_j = pipeline.b.new_job('Make toy intervals')
    in_intervals = b.read_input(str(refs.calling_interval_lists[SequencingType.EXOME]))
    intervals_j.command(f"""
    grep ^@ {in_intervals} > {intervals_j.out}
    grep -v ^@  {in_intervals} \
    | awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01) print $0 }}' \
    >> {intervals_j.out}
    """)
    out_intervals = refs.calling_interval_lists[SequencingType.TOY]
    b.write_output(intervals_j.out, str(out_intervals))

    fasta = refs.fasta_res_group(b)
    for s in samples:
        cram_j = b.new_job('Subset CRAM', s.get_job_attrs())
        cram_j.depends_on(intervals_j)
        cram_j.image(images.SAMTOOLS_PICARD_IMAGE)
        nthreads = STANDARD.set_resources(cram_j, fraction=0.5).get_nthreads()
        cram = cast(CramPath, s.alignment_input).resource_group(b)

        cram_j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            }
        )

        cram_j.command(f"""
        grep -v ^@ {intervals_j.out} > regions.bed
        
        samtools view {cram.cram} -@{nthreads - 1} \
        -L regions.bed --subsample 0.1 \
        -T {fasta.base} -Ocram -o {cram_j.output_cram.cram}
        
        samtools index -@{nthreads - 1} {cram_j.output_cram.cram} \
        {cram_j.output_cram['cram.crai']}
        """)
        b.write_output(
            cram_j.output_cram, 
            str(SUBSET_CRAM_BY_SID[s.id]).replace('.cram', '')
        )

        fastq_j = extract_fastq(
            b=b,
            cram=cram_j.output_cram,
            ext='cram',
            refs=refs,
            job_attrs=s.get_job_attrs(),
            output_fq1=SUBSET_FQ_BY_SID[s.id].r1,
            output_fq2=SUBSET_FQ_BY_SID[s.id].r2,
        )
        fastq_j.depends_on(cram_j)

    pipeline.run(wait=True)


def make_gvcfs_for_joint_calling(pipeline):
    """
    Subset GVCFs to exome for joint-calling.
    """

    b = pipeline.b
    refs = pipeline.refs
    exome_intervals = b.read_input(
        str(refs.calling_interval_lists[SequencingType.EXOME])
    )

    p = pipeline.cohort.create_dataset(DATASET)
    samples = [p.add_sample(sid, external_id=sid) for sid in SAMPLES]
    for s in samples:
        j = pipeline.b.new_job('Subset GVCF', dict(sample=s.id))
        j.image(images.BCFTOOLS_IMAGE)
        inp_gvcf = pipeline.b.read_input_group(
            **{
                'g.vcf.gz': str(FULL_GVCF_BY_SID[s.id]),
                'g.vcf.gz.tbi': str(FULL_GVCF_BY_SID[s.id]) + '.tbi',
            }
        )
        j.declare_resource_group(
            output_gvcf={
                'g.vcf.gz': '{root}-' + s.id + '.g.vcf.gz',
                'g.vcf.gz.tbi': '{root}-' + s.id + '.g.vcf.gz.tbi',
            }
        )
        j.command(
            f"""
            grep -v ^@ {exome_intervals} > regions.bed

            bcftools view -R regions.bed {inp_gvcf['g.vcf.gz']} \\
            -Oz -o {j.output_gvcf['g.vcf.gz']}
            tabix -p vcf {j.output_gvcf['g.vcf.gz']}
            """
        )
        b.write_output(
            j.output_gvcf, str(SUBSET_GVCF_BY_SID[s.id]).replace('.g.vcf.gz', '')
        )

    pipeline.run(wait=True)


if __name__ == '__main__':
    main()
