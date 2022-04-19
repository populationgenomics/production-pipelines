"""
Hail query function to export an annotated matrix table into a VCF file.
"""

import logging
import hail as hl
from gnomad.utils.vep import vep_struct_to_csq

logger = logging.getLogger(__file__)


def export_vcf(mt_path: str, out_vcf_path: str):
    """
    Export a matrix table into a VCF file.
    """
    mt = hl.read_matrix_table(mt_path)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            CSQ=vep_struct_to_csq(mt.vep),
            # clinvar
            ALLELEID=mt.clinvar.allele_id,
            CLNSIG=mt.clinvar.clinical_significance,
            CLNRECSTAT=mt.clinvar.gold_stars,
            # cadd
            cadd_phred=mt.cadd.PHRED,
            # eigen
            Eigen_phred=mt.eigen.Eigen_phred,
            # dbnsfp
            SIFT_pred=mt.dbnsfp.SIFT_pred,
            Polyphen2_HVAR_pred=mt.dbnsfp.Polyphen2_HVAR_pred,
            MutationTaster_pred=mt.dbnsfp.MutationTaster_pred,
            FATHMM_pred=mt.dbnsfp.FATHMM_pred,
            MetaSVM_pred=mt.dbnsfp.MetaSVM_pred,
            REVEL_score=mt.dbnsfp.REVEL_score,
            GERP_RS=mt.dbnsfp.GERP_RS,
            phastCons100way_vertebrate=mt.dbnsfp.phastCons100way_vertebrate,
            # g1k
            g1k_AC=mt.g1k.AC,
            g1k_AF=mt.g1k.AF,
            g1k_AN=mt.g1k.AN,
            g1k_POPMAX_AF=mt.g1k.POPMAX_AF,
            # geno2mp
            geno2mp_HPO_Count=mt.geno2mp.HPO_Count,
            # gnomad_genomes
            gnomad_genomes_AF=mt.gnomad_genomes.AF,
            gnomad_genomes_AN=mt.gnomad_genomes.AN,
            gnomad_genomes_AC=mt.gnomad_genomes.AC,
            gnomad_genomes_Hom=mt.gnomad_genomes.Hom,
            gnomad_genomes_AF_POPMAX_OR_GLOBAL=mt.gnomad_genomes.AF_POPMAX_OR_GLOBAL,
            gnomad_genomes_FAF_AF=mt.gnomad_genomes.FAF_AF,
            gnomad_genomes_Hemi=mt.gnomad_genomes.Hemi,
            # mpc
            MPC=mt.mpc.MPC,
            # primate_ai
            primate_ai_score=mt.primate_ai.score,
            # splice_ai
            splice_ai_delta_score=mt.splice_ai.delta_score,
            # topmed
            topmed_AC=mt.topmed.AC,
            topmed_AF=mt.topmed.AF,
            topmed_AN=mt.topmed.AN,
            topmed_Hom=mt.topmed.Hom,
            topmed_Het=mt.topmed.Het,
        )
    )
    hl.export_vcf(mt, out_vcf_path, tabix=True)
