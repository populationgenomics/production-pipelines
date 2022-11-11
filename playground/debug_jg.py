#!/usr/bin/env python3

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import image_path, fasta_res_group, command, reference_path
from cpg_workflows import get_batch, get_cohort
from cpg_workflows.jobs import vep
from cpg_workflows.resources import STANDARD

b = get_batch(
    'Debug exomes missing variant https://github.com/populationgenomics/seqr-private/issues/32'
)
# reference = fasta_res_group(b)
# tar_path = 'gs://cpg-seqr-main/exome/JointGenotyping/dd7b2003026c7a6c70057a9c0f170074be6322_628/genomicsdbs/interval_60_outof_100.tar'
# tar = b.read_input(tar_path)
# interval = b.read_input(
#     'gs://cpg-seqr-main/exome/JointGenotyping/dd7b2003026c7a6c70057a9c0f170074be6322_628/intervals_100/60.interval_list'
# )


def _run_genotypegvcfs():
    j = b.new_job('GenotypeGVCFs')
    j.image(image_path('gatk'))
    STANDARD.set_resources(j, ncpu=4)
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    tar -xf {tar} -C $BATCH_TMPDIR/
    WORKSPACE=gendb://$BATCH_TMPDIR/{to_path(tar_path).with_suffix('').name}

    gatk --java-options "-Xmx14000m" \
    GenotypeGVCFs \
    --verbosity DEBUG \
    -R {reference.base} \
    -O {j.output_vcf['vcf.gz']} \
    -D {reference_path('broad/dbsnp_vcf')} \
    -V $WORKSPACE \
    -L {interval} \
    --only-output-calls-starting-in-intervals \
    --merge-input-intervals \
    -G AS_StandardAnnotation
    """
    j.command(
        command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
    )
    b.write_output(j.output_vcf, 'gs://cpg-seqr-main-tmp/debug/joint-calling.vcf.gz')


def _extract_vcf():
    j = b.new_job('Extract VCF')
    j.image(image_path('gatk'))
    STANDARD.set_resources(j, ncpu=4)
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    mv {tar} interval_60_outof_100.tar 
    tar -xf interval_60_outof_100.tar
    WORKSPACE=gendb://interval_60_outof_100
    
    gatk --java-options "-Xmx14000m" \
    SelectVariants \
    -R {reference.base} \\
    -O {j.output_vcf['vcf.gz']} \\
    -V $WORKSPACE \
    -L chr12:25241000-25246000
    """
    j.command(
        command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
    )
    b.write_output(
        j.output_vcf,
        'gs://cpg-seqr-main-tmp/debug/extract_vcf_chr12:25241000-25246000',
    )


def _vep():
    cohort = get_cohort()
    scatter_count = 50
    if len(cohort.get_samples()) > 300:
        scatter_count = 100
    if len(cohort.get_samples()) > 1000:
        scatter_count = 200

    jobs = vep.add_vep_jobs(
        get_batch(),
        input_siteonly_vcf_path=to_path(
            'gs://cpg-seqr-main/exome/JointGenotyping/ce48505980adbb5c1b1a4bf865548fa7927992_861-siteonly.vcf.gz'
        ),
        out_path=to_path('gs://cpg-seqr-main-tmp/debug/vep.ht'),
        tmp_prefix=to_path('gs://cpg-seqr-main-tmp/debug'),
        scatter_count=scatter_count,
    )


# _run_genotypegvcfs()
# _add_joint_genotyper_job(
#     get_batch(),
#     genomicsdb_path=to_path(
#         'gs://cpg-seqr-main/JointGenotyping/6ae7ac744240e459f9e38f794631957c066ae4_1359/genomicsdbs/interval_6_outof_200.tar'
#     ),
#     overwrite=False,
#     number_of_samples=1000,
#     interval=b.read_input(
#         'gs://cpg-seqr-main-tmp/batch-tmp/ed12a2/Make_200_intervals_for_genome-1R6AN/6.interval_list'
#     ),
# )
# _vep()


from analysis_runner import dataproc


j = b.new_job('test')
j._preemptible = False
b.run(wait=False)
