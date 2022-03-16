"""
Docker images used in the pipelines.
"""
import logging

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:c27f80a5d68f5925614c95ef51b09dd0ecf15074-hail-0.2.85.post2-a6471564e3026200c8ccb6539210bca99ef80260'

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.3.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:4.0.0'
BIOINFO_IMAGE = f'{AR_REPO}/bioinformatics:v1-1'  # includes aligners
BWAMEM2_IMAGE = f'{AR_REPO}/bazam_bwamem2:v0'
SAMTOOLS_PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:v0.2.15'
PEDDY_IMAGE = f'{AR_REPO}/peddy:0.4.8--pyh5e36f6f_0'
VEP_IMAGE = f'{AR_REPO}/vep:105'
VERIFY_BAMID = f'{AR_REPO}/verify-bam-id:1.0.1'
