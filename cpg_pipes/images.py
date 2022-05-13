"""
Docker images used in the pipelines.
"""
import os

from cpg_utils.hail_batch import image_path

os.environ['CPG_IMAGE_REGISTRY_PREFIX'] = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_IMAGE = image_path('gatk:4.2.3.0')
BCFTOOLS_IMAGE = image_path('bcftools:1.10.2--h4f4756c_2')
SM_IMAGE = image_path('sm-api:4.0.0')
BWA_IMAGE = image_path('bwa:v0')
BWAMEM2_IMAGE = image_path('bwamem2:v0')
DRAGMAP_IMAGE = image_path('dragmap:1.2.1')
SAMTOOLS_PICARD_IMAGE = image_path('picard_samtools:v0')
SOMALIER_IMAGE = image_path('somalier:v0.2.15')
PEDDY_IMAGE = image_path('peddy:v0')
VEP_IMAGE = image_path('vep:105')
VERIFY_BAMID = image_path('verify-bam-id:1.0.1')
MULTIQC_IMAGE = image_path('multiqc:v1.12')
FASTQC_IMAGE = image_path('fastqc:v0.11.9_cv8')
DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:7f41eeb90c8bec8836a1cd20ad1911b5989a5893-hail-db65c33c29100c64405c39ebace90a7c463b4bec'
