"""
Reference data and indices on used in the pipelines.
"""
from os.path import join

REF_BUCKET = 'gs://cpg-reference'

# VEP
VEP_LOFTEE = join(REF_BUCKET, 'vep/loftee_GRCh38.tar')
VEP_CACHE = join(REF_BUCKET, 'vep/homo_sapiens_vep_105_GRCh38.tar')
