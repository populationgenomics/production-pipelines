#!/usr/bin/env python3

"""
Generate data for unit tests
"""
from typing import cast
from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.hail_batch import reference_path, image_path, fasta_res_group
from cpg_utils.flows.resources import STANDARD
from cpg_utils.flows.workflow import Workflow
from cpg_utils.flows.filetypes import CramPath
from jobs.align import extract_fastq
import utils


def main():
    """
    Generate data for unit tests
    """
    pipeline = Pipeline(
        name='make_test_data',
        description='Make test data',
    )
    # make_subset_crams(pipeline)
    # make_gvcfs_for_joint_calling(pipeline)
    # jointcalling_vcf_to_exome(pipeline)
    jointcalling_vcf_to_sub_exome(pipeline, fraction=0.001)


def make_subset_crams(pipeline: Pipeline):
    """
    Make toy CRAMs that span entire genome (not just chr20), but are small:
    1. Take WGS CRAMs
    2. Randomly select 1% of exome regions
    3. Reduce the coverage to 1% of original, write CRAMs
    4. Extract FASTQ pairs

    For alignment testing.
    """
    b = pipeline.b
    d = pipeline.cohort.create_dataset(utils.DATASET)
    samples = [
        d.add_sample(
            sid,
            external_id=sid,
            alignment_input_by_seq_type={
                utils.SEQ_TYPE: CramPath(utils.FULL_CRAM_BY_SID[sid])
            },
        )
        for sid in utils.SAMPLES
    ]

    intervals_j = pipeline.b.new_job('Make toy intervals: 1% of exome')
    in_intervals = b.read_input(
        str(reference_path('broad/exome_calling_interval_lists'))
    )
    intervals_j.command(
        f"""
    grep ^@ {in_intervals} > {intervals_j.out}
    grep -v ^@  {in_intervals} \
    | awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01) print $0 }}' \
    >> {intervals_j.out}
    """
    )

    fasta = fasta_res_group(b)
    for s in samples:
        cram_j = b.new_job('Subset CRAM', s.get_job_attrs())
        cram_j.depends_on(intervals_j)
        cram_j.image(image_path('samtools'))
        nthreads = STANDARD.set_resources(cram_j, fraction=0.5).get_nthreads()
        cram = cast(
            CramPath, s.alignment_input_by_seq_type[utils.SEQ_TYPE]
        ).resource_group(b)

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
            job_attrs=s.get_job_attrs(),
            output_fq1=utils.TOY_FQ_BY_SID[s.id].r1,
            output_fq2=utils.TOY_FQ_BY_SID[s.id].r2,
        )
        fastq_j.depends_on(cram_j)

    pipeline.run(wait=True)


def make_gvcfs_for_joint_calling(pipeline, pct=1):
    """
    Subset GVCFs to {pct}% of exome for joint-calling test.
    """

    b = pipeline.b
    refs = pipeline.refs

    intervals_j = pipeline.b.new_job(f'Make toy intervals: {pct}% of exome')
    in_intervals = b.read_input(str(refs.calling_interval_lists['exome']))
    intervals_j.command(
        f"""
    grep ^@ {in_intervals} > {intervals_j.out}
    grep -v ^@  {in_intervals} \
    | awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= {pct / 100}) print $0 }}' \
    >> {intervals_j.out}
    """
    )
    out_intervals_path = (
        utils.BASE_BUCKET / f'inputs/exome{pct}pct/calling_regions.interval_list'
    )
    b.write_output(intervals_j.out, str(out_intervals_path))

    d = pipeline.cohort.create_dataset(utils.DATASET)
    samples = [d.add_sample(sid, external_id=sid) for sid in utils.SAMPLES]
    for s in samples:
        _subset_vcf(
            b,
            utils.FULL_GVCF_BY_SID[s.run_id],
            out_intervals_path,
            utils.EXOME_1PCT_GVCF_BY_SID[s.run_id],
        )

    pipeline.run(wait=True)


def jointcalling_vcf_to_exome(pipeline):
    """
    Subset joint-calling VCF to exome.
    """
    b = pipeline.b
    refs = pipeline.refs

    # Subsetting full to exomes
    _subset_vcf(
        b,
        utils.BASE_BUCKET / 'inputs/full/9samples-joint-called.vcf.gz',
        refs.calling_interval_lists['exome'],
        utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called.vcf.gz',
    )
    _subset_vcf(
        b,
        utils.BASE_BUCKET / 'inputs/full/9samples-joint-called-siteonly.vcf.gz',
        refs.calling_interval_lists['exome'],
        utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called-siteonly.vcf.gz',
    )
    pipeline.run(wait=True)


def jointcalling_vcf_to_sub_exome(pipeline, fraction=0.05):
    """
    Subset joint-calling VCF to a fraction of exome.
    """
    b = pipeline.b

    intervals_j = pipeline.b.new_job(f'Make toy intervals: {fraction} of exome')
    in_intervals = b.read_input(
        str(reference_path('broad/exome_calling_interval_lists'))
    )
    intervals_j.command(
        f"""
    grep ^@ {in_intervals} > {intervals_j.out}
    grep -v ^@  {in_intervals} \
    | awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= {fraction}) print $0 }}' \
    >> {intervals_j.out}
    """
    )
    pct = fraction * 100
    out_intervals_path = (
        utils.BASE_BUCKET / f'inputs/exome{pct}pct/calling_regions.interval_list'
    )
    b.write_output(intervals_j.out, str(out_intervals_path))

    _subset_vcf(
        b,
        utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called.vcf.gz',
        out_intervals_path,
        utils.BASE_BUCKET / f'inputs/exome{pct}pct/9samples-joint-called.vcf.gz',
    )
    _subset_vcf(
        b,
        utils.BASE_BUCKET / 'inputs/exome/9samples-joint-called-siteonly.vcf.gz',
        out_intervals_path,
        utils.BASE_BUCKET
        / f'inputs/exome{pct}pct/9samples-joint-called-siteonly.vcf.gz',
    )
    pipeline.run(wait=True)


def _subset_vcf(
    b: Batch,
    vcf_path: Path,
    intervals_path: Path,
    out_path: Path,
) -> Job:
    """
    Make job that subsets VCF or GVCF to intervals.
    """
    intervals = b.read_input(str(intervals_path))

    j = b.new_job(f'Subset VCF {vcf_path} -> {out_path}')
    j.image(image_path('bcftools'))

    vcf = b.read_input_group(
        **{
            'vcf.gz': str(vcf_path),
            'vcf.gz.tbi': str(vcf_path) + '.tbi',
        }
    )
    j.declare_resource_group(
        out_vcf={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        }
    )
    j.command(
        f"""
        grep -v ^@ {intervals} > regions.bed

        bcftools view -R regions.bed {vcf['vcf.gz']} \\
        -Oz -o {j.out_vcf['vcf.gz']}
        tabix -p vcf {j.out_vcf['vcf.gz']}
        """
    )
    STANDARD.set_resources(j, fraction=0.5)
    b.write_output(j.out_vcf, str(out_path).replace('.vcf.gz', ''))
    return j


if __name__ == '__main__':
    main()
