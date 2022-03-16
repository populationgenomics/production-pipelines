"""
Reference data and indices on used in the pipelines.
"""
from os.path import join
import logging

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

REF_BUCKET = 'gs://cpg-reference'
BROAD_REF_BUCKET = f'{REF_BUCKET}/hg38/v1'

# BED files
NOALT_REGIONS = join(REF_BUCKET, 'noalt.bed')

# Somalier
SOMALIER_SITES = join(REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')
SOMALIER_1KG_TARGZ = join(REF_BUCKET, 'somalier/v0/1kg.somalier.tar.gz')
SOMALIER_1KG_LABELS_TSV = join(REF_BUCKET, 'somalier/v0/ancestry-labels-1kg.tsv')

# VEP
VEP_LOFTEE = join(REF_BUCKET, 'vep/loftee_GRCh38.tar')
VEP_CACHE = join(REF_BUCKET, 'vep/homo_sapiens_vep_105_GRCh38.tar')

# Fasta
REF_FASTA = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.fasta')
REF_D = dict(
    base=REF_FASTA,
    fai=REF_FASTA + '.fai',
    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
    + '.dict',
)


# DRAGMAP indices
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

# BWA indices
BWA_INDEX_EXTS = ['sa', 'amb', 'bwt', 'ann', 'pac', 'alt']
BWAMEM2_INDEX_EXTS = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac', 'alt']
BWAMEM2_INDEX_PREFIX = REF_FASTA

# GATK intervals
UNPADDED_INTERVALS = join(BROAD_REF_BUCKET, 'hg38.even.handcurated.20k.intervals')
NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS = 50
NUMBER_OF_GENOMICS_DB_INTERVALS = 50
PRECOMPUTED_INTERVALS = {
    50: join(REF_BUCKET, 'intervals', '50intervals'),
    10: join(REF_BUCKET, 'intervals', '10intervals'),
}

# VQSR
DBSNP_VCF = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf')
DBSNP_VCF_INDEX = join(BROAD_REF_BUCKET, 'Homo_sapiens_assembly38.dbsnp138.vcf.idx')
HAPMAP_RESOURCE_VCF = join(BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz')
HAPMAP_RESOURCE_VCF_INDEX = join(BROAD_REF_BUCKET, 'hapmap_3.3.hg38.vcf.gz.tbi')
OMNI_RESOURCE_VCF = join(BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz')
OMNI_RESOURCE_VCF_INDEX = join(BROAD_REF_BUCKET, '1000G_omni2.5.hg38.vcf.gz.tbi')
ONE_THOUSAND_GENOMES_RESOURCE_VCF = join(BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz')
ONE_THOUSAND_GENOMES_RESOURCE_VCF_INDEX = join(BROAD_REF_BUCKET, '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi')
MILLS_RESOURCE_VCF = join(BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz')
MILLS_RESOURCE_VCF_INDEX = join(BROAD_REF_BUCKET, 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi')
AXIOM_POLY_RESOURCE_VCF = join(BROAD_REF_BUCKET, 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz')
AXIOM_POLY_RESOURCE_VCF_INDEX = join(BROAD_REF_BUCKET, 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi')

# Contamination check
CONTAM_BUCKET = f'{BROAD_REF_BUCKET}/contamination-resources/1000g'
WGS_COVERAGE_INTERVAL_LIST = f'{BROAD_REF_BUCKET}/wgs_coverage_regions.hg38.interval_list'
