"""
Docker images used in the pipelines.
"""

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:7f41eeb90c8bec8836a1cd20ad1911b5989a5893-hail-db65c33c29100c64405c39ebace90a7c463b4bec'

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.3.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:4.0.0'
BWA_IMAGE = f'{AR_REPO}/bwa:v0'
BWAMEM2_IMAGE = f'{AR_REPO}/bwamem2:v0'
DRAGMAP_IMAGE = f'{AR_REPO}/dragmap:1.2.1'
SAMTOOLS_PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:v0.2.15'
PEDDY_IMAGE = f'{AR_REPO}/peddy:v0'
VEP_IMAGE = f'{AR_REPO}/vep:105'
VERIFY_BAMID = f'{AR_REPO}/verify-bam-id:1.0.1'
MULTIQC_IMAGE = f'{AR_REPO}/multiqc:v1.12'
FASTQC_IMAGE = f'{AR_REPO}/fastqc:v0.11.9_cv8'

# For exploration
BIOINFO_IMAGE = f'{AR_REPO}/bioinformatics:v1-1'
