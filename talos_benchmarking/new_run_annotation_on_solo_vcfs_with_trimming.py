"""
takes the ms VCFs we have available, and runs the Talos annotation pipeline on them.

This version takes the single-sample VCFs, but passes each through a 'variants only in this sample' filter
"""

import logging
import random

from cpg_utils import hail_batch, to_path


logging.basicConfig(level=logging.INFO)


random.seed(42)

batch_instance = hail_batch.get_batch('SingleSample Merge Talos Nextflow')

# read some reference files to feed into the annotation pipeline
echtvar_gnomad = "gs://cpg-common-main/gnomad/echtvar/gnomad_4.1_region_merged_GRCh38_whole_genome.zip"
echtvar_local = hail_batch.get_batch('Multisample Talos Nextflow').read_input(echtvar_gnomad)

mane_json = "gs://cpg-common-main/references/mane_1.4/mane_1.4.json"
mane_local = hail_batch.get_batch().read_input(mane_json)

ensembl_bed = "gs://cpg-common-main/references/ensembl_113/GRCh38.bed"
ensembl_local = hail_batch.get_batch().read_input(ensembl_bed)

merged_bed = 'gs://cpg-common-main/references/ensembl_113/merged_GRCh38.bed'
merged_bed_local = hail_batch.get_batch().read_input(merged_bed)

unmasked_reference = "gs://cpg-common-main/references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta"
reference_local = hail_batch.get_batch().read_input(unmasked_reference)

am_tar = 'gs://cpg-common-test/references/alphamissense/alphamissense_38.ht.tar'
am_local = hail_batch.get_batch().read_input(am_tar)

ensembl_gff = 'gs://cpg-common-main/references/ensembl_113/GRCh38.gff3.gz'
gff_local = hail_batch.get_batch().read_input(ensembl_gff)

# where are the solo VCFs? This is a combination of GHFM-kidgen and acute-care
GS_VCFS = 'gs://cpg-seqr-test/talos_benchmarking/solo_vcfs'

# grab all the VCFs in the solo_vcfs bucket
vcf_list = [str(each_vcf) for each_vcf in to_path(GS_VCFS).glob('*.vcf.gz')]

vcf_inputs = [batch_instance.read_input_group(gvcf=sample_vcf, index=f'{sample_vcf}.tbi') for sample_vcf in vcf_list]

image = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/talos:PR_552'

for each_count in [5, 10, 25, 50, 100, 250, 375, 600, 1000]:

    logging.info(f'Considering batch size {each_count}')

    output_folder = f'gs://cpg-acute-care-test/talos_benchmarking/new_trimmed_ms_merged_results/{each_count}'

    if to_path(f'{output_folder}/report.html').exists():
        logging.info(f'{output_folder}/report.html exists, skipping')
        continue

    new_job = hail_batch.get_batch().new_bash_job(f'Run Nextflow for {each_count} MS VCF')

    new_job.cpu(16).memory('32GiB').storage(f'{min(each_count, 500)}GiB')
    new_job.image(image)

    # create a subset of VCFs to run
    vcf_group = random.sample(vcf_inputs, min(each_count, len(vcf_inputs)))

    logging.info(f'Detected {len(vcf_group)} input vcfs, of the expected {each_count}')

    new_job.command('set -ex')
    new_job.command('echo "Inputs copied in"')
    new_job.command('date')
    new_job.command('mkdir $BATCH_TMPDIR/individual_vcfs')
    new_job.command('mkdir $BATCH_TMPDIR/work')

    # move these into --input_vcf_dir, via a bcftools trim (might take way longer?)
    for enum_count, each_vcf in enumerate(vcf_group):
        new_job.command(f'bcftools view -c1 -W=tbi -Oz -o $BATCH_TMPDIR/individual_vcfs/{enum_count}.vcf.gz {each_vcf.gvcf}')
        new_job.command(f'rm {each_vcf.gvcf} {each_vcf.index}')

    new_job.command('echo "VCFs finished processing"')
    new_job.command('date')

    new_job.command(
        f"""
    mkdir $BATCH_TMPDIR/output

    nextflow -log {new_job.log} -c nextflow/annotation.config run nextflow/annotation.nf \\
        -without-docker -with-report {new_job.report} -w $BATCH_TMPDIR/work \\
        --input_vcf_dir $BATCH_TMPDIR/individual_vcfs \\
        --alphamissense_tar {am_local} \\
        --cohort {each_count} \\
        --cohort_output_dir $BATCH_TMPDIR/output \\
        --ensembl_bed {ensembl_local} \\
        --ensembl_merged_bed {merged_bed_local} \\
        --ensembl_gff {gff_local} \\
        --gnomad_zip {echtvar_local} \\
        --mane_json {mane_local} \\
        --ref_genome {reference_local}

    gcloud storage cp -r $BATCH_TMPDIR/output/{each_count}.mt {output_folder}/
    """
    )
    new_job.command('echo "Completed"')
    new_job.command('date')

    hail_batch.get_batch().write_output(new_job.log, f'{output_folder}/nextflow.log')
    hail_batch.get_batch().write_output(new_job.report, f'{output_folder}/report.html')

hail_batch.get_batch().run(wait=False)
