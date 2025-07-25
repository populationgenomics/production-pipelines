"""
takes the ms VCFs we have available, and runs the Talos annotation pipeline on them.
"""

import logging

from cpg_utils import hail_batch, to_path


logging.basicConfig(level=logging.INFO)


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

for each_group in [5, 10, 25, 50, 100, 250]:

    output_folder = f'gs://cpg-acute-care-test/talos_benchmarking/ms_results_trimmed/{each_group}'

    if to_path(f'{output_folder}/report.html').exists():
        logging.info(f'{output_folder}/report.html exists, skipping')
        continue

    ms_vcf = f'gs://cpg-acute-care-test/talos_benchmarking/ms_vcfs/{each_group}.vcf.bgz'
    if not to_path(ms_vcf).exists():
        logging.info(f'skipping {each_group}, input doesn\'t exist yet')
        continue

    local_vcf = hail_batch.get_batch().read_input_group(gvcf=ms_vcf, index=f'{ms_vcf}.tbi').gvcf

    image = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/talos:PR_552'

    new_job = hail_batch.get_batch().new_bash_job(f'Run Nextflow for {each_group} MS VCF')
    new_job.cpu(16).memory('32GiB').storage('250GiB')
    new_job.image(image)

    new_job.command(f"""
    set -x
    
    mkdir $BATCH_TMPDIR/output
    
    nextflow -log {new_job.log} -c nextflow/annotation.config run nextflow/annotation.nf \\
        -without-docker -with-report {new_job.report} \\
        --merged_vcf {local_vcf} \\
        --alphamissense_tar {am_local} \\
        --cohort {each_group} \\
        --cohort_output_dir $BATCH_TMPDIR/output \\
        --ensembl_bed {ensembl_local} \\
        --ensembl_merged_bed {merged_bed_local} \\
        --ensembl_gff {gff_local} \\
        --gnomad_zip {echtvar_local} \\
        --mane_json {mane_local} \\
        --ref_genome {reference_local}
    
    gcloud storage cp -r $BATCH_TMPDIR/output/{each_group}.mt {output_folder}/
    """)

    hail_batch.get_batch().write_output(new_job.log, f'{output_folder}/nextflow.log')
    hail_batch.get_batch().write_output(new_job.report, f'{output_folder}/report.html')

hail_batch.get_batch().run(wait=False)
