"""
Utility module to encapsulate pointers to the reference files used in the pipeline
"""

from os.path import join


# Images
DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:ca793d8df6dd8fa20f1908d5c76eb98ab2cf53e4-hail-0.2.73.devc6f6f09cec08'

AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
GATK_VERSION = '4.2.3.0'
GATK_IMAGE = f'{AR_REPO}/gatk:{GATK_VERSION}'
BCFTOOLS_IMAGE = f'{AR_REPO}/bcftools:1.10.2--h4f4756c_2'
SM_IMAGE = f'{AR_REPO}/sm-api:4.0.0'
BIOINFO_IMAGE = f'{AR_REPO}/bioinformatics:v1-0'  # include aligners
BWAMEM2_IMAGE = f'{AR_REPO}/bazam_bwamem2:v0'
SAMTOOLS_PICARD_IMAGE = f'{AR_REPO}/picard_samtools:v0'
SOMALIER_IMAGE = f'{AR_REPO}/somalier:latest'
PEDDY_IMAGE = f'{AR_REPO}/peddy:0.4.8--pyh5e36f6f_0'


# Files
REF_BUCKET = 'gs://cpg-reference'
NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')
SOMALIER_SITES = join(REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')

BROAD_REF_BUCKET = f'{REF_BUCKET}/hg38/v1'

REF_FASTA = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
REF_D = dict(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
)

DRAGMAP_INDEX_BUCKET = f'{REF_BUCKET}/dragmap/v0'
DRAGMAP_INDEX_FILES = [
    'hash_table.cfg',
    'hash_table.cfg.bin',
    'hash_table.cmp',
    'hash_table_stats.txt',
    'ref_index.bin',
    'reference.bin',
    'repeat_mask.bin',
    'str_table.bin',
]

BWAMEM2_INDEX_PREFIX = REF_FASTA
BWAMEM2_INDEX_EXTS = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac']
BWA_INDEX_EXTS = ['sa', 'amb', 'bwt', 'ann', 'pac']

DBSNP_VCF = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
UNPADDED_INTERVALS = join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')

NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50
PRECOMPUTED_INTERVALS = {
    50: join(REF_BUCKET, 'intervals', '50intervals'),
    10: join(REF_BUCKET, 'intervals', '10intervals'),
}
