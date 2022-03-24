#!/usr/bin/env python3

"""
Generate data for unit tests
"""

from cpg_pipes import images, Namespace
from cpg_pipes.pipeline.pipeline import Pipeline
from .utils import DATASET, SAMPLES, FULL_GVCF_BY_SID, SUBSET_GVCF_BY_SID, setup_env


def make_gvcfs():
    """
    Subset GVCFs to chr20 for joint-calling
    """
    setup_env()
    
    pipeline = Pipeline(
        name='make_test_data',
        description='Make test data',
        analysis_dataset=DATASET,
        version='v0',
        namespace=Namespace.TEST,
    )

    jobs = []
    
    p = pipeline.cohort.add_dataset(DATASET)
    samples = [
        p.add_sample(sid, external_id=sid) for sid in SAMPLES
    ]
    for s in samples:
        subset_j = pipeline.b.new_job(
            'Subset GVCF', dict(sample_name=s.id)
        )
        subset_j.image(images.BCFTOOLS_IMAGE)
        inp_gvcf = pipeline.b.read_input_group(**{
            'g.vcf.gz': FULL_GVCF_BY_SID[s.id],
            'g.vcf.gz.tbi': FULL_GVCF_BY_SID[s.id] + '.tbi',
        })
        subset_j.declare_resource_group(
            output_gvcf={
                'g.vcf.gz': '{root}-' + s.id + '.g.vcf.gz',
                'g.vcf.gz.tbi': '{root}-' + s.id + '.g.vcf.gz.tbi',
            }
        )
        subset_j.command(f"""
            bcftools view -r chr20 {inp_gvcf['g.vcf.gz']} \\
            -Oz -o {subset_j.output_gvcf['g.vcf.gz']}
            tabix -p vcf {subset_j.output_gvcf['g.vcf.gz']}
        """)
        pipeline.b.write_output(
            subset_j.output_gvcf,
            SUBSET_GVCF_BY_SID[s.id].replace('.g.vcf.gz', '')
        )
        jobs.append(subset_j)
    
    pipeline.submit_batch(wait=True)     


make_gvcfs()
