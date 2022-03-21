"""
Reference files and indices used in bionformatics pipelines.
"""

from . import Path


class RefData:
    """
    Bioinformatics reference files, indices, and constants.
    """
    number_of_haplotype_caller_intervals = 50
    number_of_genomics_db_intervals = 50

    def fasta_res_group(self, b, indices: list | None = None):
        """
        Hail Batch resource group for fasta reference files.
        @param b: Hail Batch object.
        @param indices: list of extentions to add to the base fasta file path.
        """
        d = dict(
            base=str(self.ref_fasta),
            fai=str(self.ref_fasta) + '.fai',
            dict=str(self.ref_fasta.with_suffix('.dict')),
        )
        if indices:
            d |= {ext: f'{self.ref_fasta}.{ext}' for ext in indices}
        return b.read_input_group(**d)

    def __init__(self, bucket: Path):
        self.bucket = bucket
        
        # BED files
        self.noalt_regions = self.bucket / 'noalt.bed'
        
        # Somalier
        self.somalier_sites = self.bucket / 'somalier/v0/sites.hg38.vcf.gz'
        self.somalier_1kg_targz = self.bucket / 'somalier/v0/1kg.somalier.tar.gz'
        self.somalier_1kg_labels_tsv = self.bucket / 'somalier/v0/ancestry-labels-1kg.tsv'
    
        # VEP
        self.vep_loftee = self.bucket / 'vep/loftee_GRCh38.tar'
        self.vep_cache = self.bucket / 'vep/homo_sapiens_vep_105_GRCh38.tar' 

        self.broad_ref_bucket = bucket / 'hg38' / 'v1'

        # Fasta
        self.ref_fasta = self.broad_ref_bucket / 'Homo_sapiens_assembly38.fasta'
    
        # DRAGMAP indices
        self.dragmap_index_bucket = self.bucket / 'dragmap/v0'
        self.dragmap_index_files = [
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
        self.bwa_index_exts = ['sa', 'amb', 'bwt', 'ann', 'pac', 'alt']
        self.bwamem2_index_exts = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac', 'alt']
        self.bwamem2_index_prefix = self.ref_fasta
    
        # GATK intervals
        self.unpadded_intervals = (
            self.broad_ref_bucket / 'hg38.even.handcurated.20k.intervals'
        )
        self.precomputed_intervals = {
            50: self.bucket / 'intervals' / '50intervals',
            10: self.bucket / 'intervals' / '10intervals',
        }
    
        # VQSR
        self.dbsnp_vcf = self.broad_ref_bucket / 'Homo_sapiens_assembly38.dbsnp138.vcf'
        self.dbsnp_vcf_index = self.broad_ref_bucket / 'Homo_sapiens_assembly38.dbsnp138.vcf.idx'
        self.hapmap_resource_vcf = self.broad_ref_bucket / 'hapmap_3.3.hg38.vcf.gz'
        self.hapmap_resource_vcf_index = self.broad_ref_bucket / 'hapmap_3.3.hg38.vcf.gz.tbi'
        self.omni_resource_vcf = self.broad_ref_bucket / '1000G_omni2.5.hg38.vcf.gz'
        self.omni_resource_vcf_index = self.broad_ref_bucket / '1000G_omni2.5.hg38.vcf.gz.tbi'
        self.one_thousand_genomes_resource_vcf = self.broad_ref_bucket / '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
        self.one_thousand_genomes_resource_vcf_index = self.broad_ref_bucket / '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
        self.mills_resource_vcf = self.broad_ref_bucket / 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
        self.mills_resource_vcf_index = self.broad_ref_bucket / 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
        self.axiom_poly_resource_vcf = self.broad_ref_bucket / 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
        self.axiom_poly_resource_vcf_index = self.broad_ref_bucket / 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi'
    
        # Contamination check
        self.contam_bucket = self.broad_ref_bucket / 'contamination-resources/1000g'
        self.wgs_coverage_interval_list = self.broad_ref_bucket / 'wgs_coverage_regions.hg38.interval_list'
    
        ######################################################
        # Hail Query #
        ######################################################
        self.gnomad_ref_bucket = self.bucket / 'gnomad/v0'
        self.tel_and_cent_ht = self.gnomad_ref_bucket / 'telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht'
        self.lcr_intervals_ht = self.gnomad_ref_bucket / 'lcr_intervals/LCRFromHengHg38.ht'
        self.seg_dup_intervals_ht = self.gnomad_ref_bucket / 'seg_dup_intervals/GRCh38_segdups.ht'
        self.clinvar_ht = self.gnomad_ref_bucket / 'clinvar/clinvar_20190923.ht'
        self.hapmap_ht = self.gnomad_ref_bucket / 'hapmap/hapmap_3.3.hg38.ht'
        self.kgp_omni_ht = self.gnomad_ref_bucket / 'kgp/1000G_omni2.5.hg38.ht'
        self.kgp_hc_ht = self.gnomad_ref_bucket / 'kgp/1000G_phase1.snps.high_confidence.hg38.ht'
        self.mills_ht = self.gnomad_ref_bucket / 'mills/Mills_and_1000G_gold_standard.indels.hg38.ht'
