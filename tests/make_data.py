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

import utils


def main():
    """
    Generate data for unit tests
    """
    utils.setup_env()
    pipeline = create_pipeline(
        name='make_test_data',
        description='Make test data',
        analysis_dataset=utils.DATASET,
        namespace=Namespace.TEST,
    )
    # make_subset_crams(pipeline)
    # make_gvcfs_for_joint_calling(pipeline)
    make_joint_calling_vcf(pipeline)


def make_subset_crams(pipeline: Pipeline):
    """
    Make toy CRAMs that span entire genome, not just chr20.
    1. Take WGS CRAMs
    2. Randomly select 1% of exome regions
    3. Reduce the coverage to 1% of original, write CRAMs
    4. Extract FASTQ pairs
    """
    b = pipeline.b
    d = pipeline.cohort.create_dataset(utils.DATASET)
    samples = [
        d.add_sample(
            sid, external_id=sid, alignment_input=CramPath(utils.FULL_CRAM_BY_SID[sid])
        )
        for sid in utils.SAMPLES
    ]

    refs = pipeline.refs

    intervals_j = pipeline.b.new_job('Make toy intervals')
    in_intervals = b.read_input(str(refs.calling_interval_lists[SequencingType.EXOME]))
    intervals_j.command(
        f"""
    grep ^@ {in_intervals} > {intervals_j.out}
    grep -v ^@  {in_intervals} \
    | awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01) print $0 }}' \
    >> {intervals_j.out}
    """
    )
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

        cram_j.command(
            f"""
        grep -v ^@ {intervals_j.out} > regions.bed
        
        samtools view {cram.cram} -@{nthreads - 1} \
        -L regions.bed --subsample 0.1 \
        -T {fasta.base} -Ocram -o {cram_j.output_cram.cram}
        
        samtools index -@{nthreads - 1} {cram_j.output_cram.cram} \
        {cram_j.output_cram['cram.crai']}
        """
        )
        b.write_output(
            cram_j.output_cram, str(utils.TOY_CRAM_BY_SID[s.id]).replace('.cram', '')
        )

        fastq_j = extract_fastq(
            b=b,
            cram=cram_j.output_cram,
            ext='cram',
            refs=refs,
            job_attrs=s.get_job_attrs(),
            output_fq1=utils.TOY_FQ_BY_SID[s.id].r1,
            output_fq2=utils.TOY_FQ_BY_SID[s.id].r2,
        )
        fastq_j.depends_on(cram_j)

    pipeline.run(wait=True)


def make_gvcfs_for_joint_calling(pipeline):
    """
    Subset GVCFs to exome for joint-calling.
    UPD: not needed, replaced with make_joint_calling_vcf
    (jointgenotying works fine, it's VQSR that needs mocking)
    """

    b = pipeline.b
    refs = pipeline.refs
    exome_intervals = b.read_input(
        str(refs.calling_interval_lists[SequencingType.EXOME])
    )

    d = pipeline.cohort.create_dataset(utils.DATASET)
    samples = [d.add_sample(sid, external_id=sid) for sid in utils.SAMPLES]
    for s in samples:
        j = pipeline.b.new_job('Subset GVCF', dict(sample=s.id))
        j.image(images.BCFTOOLS_IMAGE)
        inp_gvcf = pipeline.b.read_input_group(
            **{
                'g.vcf.gz': str(utils.FULL_GVCF_BY_SID[s.id]),
                'g.vcf.gz.tbi': str(utils.FULL_GVCF_BY_SID[s.id]) + '.tbi',
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
            j.output_gvcf, str(utils.EXOME_GVCF_BY_SID[s.id]).replace('.g.vcf.gz', '')
        )

    pipeline.run(wait=True)


def make_joint_calling_vcf(pipeline):
    """
    Subset joint-calling VCF to exome.
    """
    b = pipeline.b
    refs = pipeline.refs
    exome_intervals = b.read_input(
        str(refs.calling_interval_lists[SequencingType.EXOME])
    )

    j = b.new_job('Subset joint-calling VCF')
    j.image(images.BCFTOOLS_IMAGE)

    full_vcf = utils.BASE_BUCKET / 'inputs/full/9samples-joint-called.vcf.gz'
    full_siteonly_vcf = (
        utils.BASE_BUCKET / 'inputs/full/9samples-joint-called-siteonly.vcf.gz'
    )

    vcf = b.read_input_group(
        **{
            'vcf.gz': str(full_vcf),
            'vcf.gz.tbi': str(full_vcf) + '.tbi',
        }
    )
    siteonly_vcf = b.read_input_group(
        **{
            'vcf.gz': str(full_siteonly_vcf),
            'vcf.gz.tbi': str(full_siteonly_vcf) + '.tbi',
        }
    )
    j.declare_resource_group(
        out_vcf={
            'vcf.gz': '{root}-joint-called.vcf.gz',
            'vcf.gz.tbi': '{root}-joint-called.vcf.gz.tbi',
        }
    )
    j.declare_resource_group(
        out_siteonly_vcf={
            'vcf.gz': '{root}-joint-called-siteonly.vcf.gz',
            'vcf.gz.tbi': '{root}-joint-called-siteonly.vcf.gz.tbi',
        }
    )
    j.command(
        f"""
        grep -v ^@ {exome_intervals} > regions.bed
    
        bcftools view -R regions.bed {vcf['vcf.gz']} \\
        -Oz -o {j.out_vcf['vcf.gz']}
        tabix -p vcf {j.out_vcf['vcf.gz']}

        bcftools view -R regions.bed {siteonly_vcf['vcf.gz']} \\
        -Oz -o {j.out_siteonly_vcf['vcf.gz']}
        tabix -p vcf {j.out_siteonly_vcf['vcf.gz']}
        """
    )
    STANDARD.set_resources(j, fraction=0.5)
    full_vcf = utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called.vcf.gz'
    full_siteonly_vcf = (
        utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called-siteonly.vcf.gz'
    )

    b.write_output(j.out_vcf, str(full_vcf).replace('.vcf.gz', ''))
    b.write_output(j.out_siteonly_vcf, str(full_siteonly_vcf).replace('.vcf.gz', ''))
    pipeline.run(wait=True)


if __name__ == '__main__':
    main()
