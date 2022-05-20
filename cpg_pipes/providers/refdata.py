"""
Reference files and indices used in bionformatics pipelines.
"""

from cpg_pipes.providers.storage import Path
from cpg_pipes.types import SequencingType


class RefData:
    """
    Bioinformatics reference files, indices, and constants.
    """

    number_of_shards_for_realignment = 10
    number_of_haplotype_caller_intervals = 50
    number_of_joint_calling_intervals = 50
    number_of_vep_intervals = 50
    genome_build = 'GRCh38'
    
    @property
    def broad_ref_prefix(self) -> Path:
        """
        Path to the Broad reference bucket,
        matching gs://gcp-public-data--broad-references/hg38/v0
        """
        return self.prefix / 'hg38/v0'

    @property
    def somalier_sites(self) -> Path:
        """
        Site list for somalier, https://github.com/brentp/somalier/releases/tag/v0.2.15
        """
        return self.prefix / 'somalier/v0/sites.hg38.vcf.gz'

    @property
    def somalier_1kg_targz(self) -> Path:
        """
        Somalier 1kg data for the "ancestry" command.
        """
        return self.prefix / 'somalier/v0/1kg.somalier.tar.gz'

    @property
    def somalier_1kg_labels(self) -> Path:
        """
        Somalier list of 1kg samples for the "ancestry" command.
        """
        return self.prefix / 'somalier/v0/ancestry-labels-1kg.tsv'

    @property
    def vep_loftee(self) -> Path:
        """
        LofTee reference tarball for VEP.
        """
        return self.prefix / 'vep' / 'loftee_GRCh38.tar'

    @property
    def vep_cache(self) -> Path:
        """
        Reference tarball for VEP.
        """
        return self.prefix / 'vep' / 'homo_sapiens_vep_105_GRCh38.tar'

    @property
    def vep_mount(self) -> Path:
        """
        Contains uncompressed VEP tarballs for mounting with cloudfuse.
        """
        return self.prefix / 'vep' / 'GRCh38'
    
    @property
    def intervals_prefix(self) -> Path:
        """
        To cache intervals.
        """
        return self.prefix / 'intervals'

    @property
    def gnomad_prefix(self) -> Path:
        """
        Prefix for reference files for gnomAD QC.
        """
        return self.prefix / 'gnomad/v0'

    def __init__(self, path_prefix: Path):
        self.prefix = path_prefix

        self.dragmap_ref_bucket = self.broad_ref_prefix / 'dragen_reference'
        self.ref_fasta = self.dragmap_ref_bucket / 'Homo_sapiens_assembly38_masked.fasta'
        self.dragmap_index_files = [
            'hash_table.cfg.bin',
            'hash_table.cmp',
            'reference.bin',
        ]

        # BWA indices
        self.bwa_index_exts = ['sa', 'amb', 'bwt', 'ann', 'pac', 'alt']
        self.bwamem2_index_exts = ['0123', 'amb', 'bwt.2bit.64', 'ann', 'pac', 'alt']
        self.bwamem2_index_prefix = self.ref_fasta

        # To clean up WGS VCFs
        self.noalt_regions = (
            self.broad_ref_prefix / 
            'sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz'
        )

        # Intervals
        self.calling_interval_lists = {
            SequencingType.GENOME: self.broad_ref_prefix
            / 'wgs_calling_regions.hg38.interval_list',
            SequencingType.EXOME: self.broad_ref_prefix
            / 'exome_calling_regions.v1.interval_list',
        }
        self.evaluation_interval_lists = {
            SequencingType.GENOME: self.broad_ref_prefix
            / 'wgs_evaluation_regions.hg38.interval_list',
            SequencingType.EXOME: self.broad_ref_prefix
            / 'exome_evaluation_regions.v1.interval_list',
        }
        self.wgs_coverage_interval_list = (
            self.broad_ref_prefix / 'wgs_coverage_regions.hg38.interval_list'
        )
        self.unpadded_intervals = (
            self.broad_ref_prefix / 'hg38.even.handcurated.20k.intervals'
        )
        self.exome_bed = (
            self.broad_ref_prefix
            / 'Homo_sapiens_assembly38.contam.exome_calling_regions.v1.bed'
        )

        # VQSR
        self.dbsnp_vcf = self.broad_ref_prefix / 'Homo_sapiens_assembly38.dbsnp138.vcf'
        self.dbsnp_vcf_index = (
            self.broad_ref_prefix / 'Homo_sapiens_assembly38.dbsnp138.vcf.idx'
        )
        self.hapmap_resource_vcf = self.broad_ref_prefix / 'hapmap_3.3.hg38.vcf.gz'
        self.hapmap_resource_vcf_index = (
            self.broad_ref_prefix / 'hapmap_3.3.hg38.vcf.gz.tbi'
        )
        self.omni_resource_vcf = self.broad_ref_prefix / '1000G_omni2.5.hg38.vcf.gz'
        self.omni_resource_vcf_index = (
            self.broad_ref_prefix / '1000G_omni2.5.hg38.vcf.gz.tbi'
        )
        self.one_thousand_genomes_resource_vcf = (
            self.broad_ref_prefix / '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
        )
        self.one_thousand_genomes_resource_vcf_index = (
            self.broad_ref_prefix / '1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
        )
        self.mills_resource_vcf = (
            self.broad_ref_prefix / 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
        )
        self.mills_resource_vcf_index = (
            self.broad_ref_prefix
            / 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
        )
        self.axiom_poly_resource_vcf = (
            self.broad_ref_prefix
            / 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
        )
        self.axiom_poly_resource_vcf_index = (
            self.broad_ref_prefix
            / 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi'
        )

        # Contamination check
        self.contam_bucket = self.broad_ref_prefix / 'contamination-resources/1000g'
        self.cont_ref_d = dict(
            ud=self.contam_bucket / '1000g.phase3.100k.b38.vcf.gz.dat.UD',
            bed=self.contam_bucket / '1000g.phase3.100k.b38.vcf.gz.dat.bed',
            mu=self.contam_bucket / '1000g.phase3.100k.b38.vcf.gz.dat.mu',
        )

        ######################################################
        # Hail Query #
        ######################################################
        self.tel_and_cent_ht = (
            self.gnomad_prefix
            / 'telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht'
        )
        self.lcr_intervals_ht = (
            self.gnomad_prefix / 'lcr_intervals/LCRFromHengHg38.ht'
        )
        self.seg_dup_intervals_ht = (
            self.gnomad_prefix / 'seg_dup_intervals/GRCh38_segdups.ht'
        )
        self.clinvar_ht = self.gnomad_prefix / 'clinvar/clinvar_20190923.ht'
        self.hapmap_ht = self.gnomad_prefix / 'hapmap/hapmap_3.3.hg38.ht'
        self.kgp_omni_ht = self.gnomad_prefix / 'kgp/1000G_omni2.5.hg38.ht'
        self.kgp_hc_ht = (
            self.gnomad_prefix / 'kgp/1000G_phase1.snps.high_confidence.hg38.ht'
        )
        self.mills_ht = (
            self.gnomad_prefix
            / 'mills/Mills_and_1000G_gold_standard.indels.hg38.ht'
        )

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
            str=str(self.ref_fasta.with_suffix('.str')),
        )
        if indices:
            d |= {ext: f'{self.ref_fasta}.{ext}' for ext in indices}
        return b.read_input_group(**d)
