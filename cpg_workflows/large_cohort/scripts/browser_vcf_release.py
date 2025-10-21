import copy
import itertools
import logging
import pickle
from copy import deepcopy
from pprint import pprint
from typing import Dict, List, Optional, Set, Tuple

from matplotlib import category

import hail as hl

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import genome_build, output_path
from cpg_workflows.large_cohort.browser_prepare import _freq_index_key, process_score_cutoffs
from cpg_workflows.utils import can_reuse
from gnomad.utils.filtering import add_filters_expr

logging.basicConfig(
    format='%(asctime)s (%(name)s %(lineno)s): %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


LEN_COMP_GLOBAL_ROWS = {
    'freq': ['freq_meta', 'freq_index_dict', 'freq_meta_sample_count'],
    'faf': ['faf_meta', 'faf_index_dict'],
}
LEN_COMP_JOINT_GLOBAL_ROWS = {
    'joint_freq': [
        'joint_freq_meta',
        'joint_freq_index_dict',
        'joint_freq_meta_sample_count',
    ],
    'joint_faf': ['joint_faf_meta', 'joint_faf_index_dict'],
}

GROUPS = ['adj', 'raw']
SEXES = ['XX', 'XY']
SUBSETS = {'v1': ['']}
SUBSETS = {
    'exomes': deepcopy(SUBSETS['v1']),
    'genomes': deepcopy(SUBSETS['v1']),
    'joint': [''],
}
IN_SILICO_ANNOTATIONS_INFO_DICT = None
VRS_FIELDS_DICT = None
GEN_ANC_NAMES = {
    'AFR': 'African / African American / African Caribbean',
    'AMR': 'Central & South American',
    'CSA': 'Central / South Asian',
    'EAS': 'East & Southeast Asian',
    'EUR': 'European',
    'FIL': 'Filipino',
    'MID': 'Middle Eastern / North African',
    'NA': 'Unclassified',
}
FAF_GEN_ANC_GROUPS = {
    'v1': ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'FIL', 'MID', 'NA'],
}
GEN_ANC_GROUPS = {
    'genomes': [
        'AFR',
        'AMR',
        'CSA',
        'EAS',
        'EUR',
        'FIL',
        'MID',
        'NA',
    ],
    'exomes': [
        'AFR',
        'AMR',
        'CSA',
        'EAS',
        'EUR',
        'FIL',
        'MID',
        'NA',
    ],
}
GEN_ANC_GROUPS['joint'] = list(set(GEN_ANC_GROUPS['exomes']) | set(GEN_ANC_GROUPS['genomes']))
GEN_ANC_GROUPS = {d: {pop: GEN_ANC_NAMES[pop] for pop in pops} for d, pops in GEN_ANC_GROUPS.items()}  # type: ignore
FAF_GEN_ANC_GROUPS = {
    'v1': ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'FIL', 'MID', 'NA'],
}
JOINT_FILTERS_INFO_DICT = {
    'exomes_filters': {'Description': "Filters' values from the exomes dataset."},
    'genomes_filters': {'Description': "Filters' values from the genomes dataset."},
}
VCF_INFO_REORDER = [
    'AC',
    'AN',
    'AF',
    'grpmax',
    'fafmax_faf95_max',
    'fafmax_faf95_max_gen_anc',
]
INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    'AS_QUALapprox',
    'AS_VarDP',
    'AS_MQ_DP',
    'AS_RAW_MQ',
    'AS_SB_TABLE',
]

JOINT_REGION_FLAGS_INFO_DICT = {
    'outside_capture_region': {
        'Description': 'Variant falls outside of the OurDNA exome capture regions.',
    },
    'outside_calling_region': {
        'Description': ('Variant falls outside of the OurDNA exome capture regions plus 150 bp padding.'),
    },
    'not_called_in_exomes': {
        'Description': 'Variant was not called in the OurDNA exomes.',
    },
    'not_called_in_genomes': {
        'Description': 'Variant was not called in the OurDNA genomes.',
    },
}

SORT_ORDER = [
    'subset',
    'downsampling',
    'grpmax',
    'popmax',
    'gen_anc',
    'pop',
    'subgrp',
    'subpop',
    'sex',
    'group',
]

JOINT_REGION_FLAG_FIELDS = [
    'outside_capture_region',
    'outside_calling_region',
    'not_called_in_exomes',
    'not_called_in_genomes',
]
REGION_FLAG_FIELDS_FLAT = ['lcr', 'non_par', 'segdup']  # 'nonpar' and 'decoy' are not in our region_flags field
REGION_FLAG_FIELDS = {
    'exomes': REGION_FLAG_FIELDS_FLAT
    + [
        'outside_capture_region',
        'outside_calling_region',
    ],
    'genomes': REGION_FLAG_FIELDS_FLAT,
    'joint': JOINT_REGION_FLAG_FIELDS,
}

SITE_FIELDS_FLAT = [
    'FS',
    'MQ',
    'MQRankSum',
    'QUALapprox',
    'QD',
    'ReadPosRankSum',
    'SB',
    'SOR',
    'VarDP',
]
SITE_FIELDS = {
    'exomes': SITE_FIELDS_FLAT,
    'genomes': SITE_FIELDS_FLAT,
    'joint': SITE_FIELDS_FLAT,
}
AS_FIELDS = [
    'AS_FS',
    'AS_MQ',
    'AS_MQRankSum',
    'AS_pab_max',
    'AS_QUALapprox',
    'AS_QD',
    'AS_ReadPosRankSum',
    'AS_SB_TABLE',
    'AS_SOR',
    'AS_VarDP',
]
AS_VQSR_FIELDS = ['AS_culprit', 'AS_VQSLOD']

VQSR_FIELDS = AS_VQSR_FIELDS + ['NEGATIVE_TRAIN_SITE', 'POSITIVE_TRAIN_SITE']

# These fields are not in our genome/exome frequency tables - they come from vqsr_ht.allele_info
ALLELE_TYPE_FIELDS_FLAT = [
    'allele_type',
    'has_star',
    'n_alt_alleles',
    'variant_type',
    'was_mixed',
]
ALLELE_TYPE_FIELDS = {
    'exomes': ALLELE_TYPE_FIELDS_FLAT,
    'genomes': ALLELE_TYPE_FIELDS_FLAT,
    'joint': ALLELE_TYPE_FIELDS_FLAT,
}

HISTS = ['gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt']

CURRENT_VEP_VERSION = '110'
VEP_CSQ_FIELDS = {
    '101': 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info',
    '105': 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|UNIPROT_ISOFORM|SOURCE|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info',
    '110': 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|UNIPROT_ISOFORM|SOURCE|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info|EXISTING_UORFS|EXISTING_INFRAME_OORFS|EXISTING_OUTOFFRAME_OORFS|5UTR_CONSEQUENCE|5UTR_ANNOTATION|AM_CLASS|AM_PATHOGENICITY',
}

VEP_CSQ_HEADER = 'Consequence annotations from Ensembl VEP. Format:' f' {VEP_CSQ_FIELDS[CURRENT_VEP_VERSION]}'

####
# From our large cohorts vqsr.py
STANDARD_FEATURES = [
    'ReadPosRankSum',
    'MQRankSum',
    'QD',
    'FS',
    'SOR',
]
SNP_STANDARD_FEATURES = STANDARD_FEATURES + ['MQ']
INDEL_STANDARD_FEATURES = STANDARD_FEATURES

ALLELE_SPECIFIC_FEATURES = [
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    # Not using depth for the following reasons:
    # 1. The Broad pipelines don't use it;
    # 2. -G AS_StandardAnnotation flag to GenotypeGVCFs doesn't include it;
    # 3. For exomes, depth is an irrelevant feature and should be skipped:
    # 'AS_VarDP'
    # Note that for consistency, we also skip it for WGS.
]
SNP_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES + ['AS_MQ']
INDEL_ALLELE_SPECIFIC_FEATURES = ALLELE_SPECIFIC_FEATURES
####

INFO_DICT = {
    'FS': {
        'Description': 'Phred-scaled p-value of Fisher"s exact test for strand bias',
    },
    'MQ': {
        'Description': ('Root mean square of the mapping quality of reads across all samples'),
    },
    'MQRankSum': {
        'Description': ('Z-score from Wilcoxon rank sum test of alternate vs. reference read' ' mapping qualities'),
    },
    'QD': {
        'Description': ('Variant call confidence normalized by depth of sample reads supporting a' ' variant'),
    },
    'ReadPosRankSum': {
        'Description': ('Z-score from Wilcoxon rank sum test of alternate vs. reference read' ' position bias'),
    },
    'SOR': {'Description': 'Strand bias estimated by the symmetric odds ratio test'},
    'POSITIVE_TRAIN_SITE': {
        'Description': ('Variant was used to build the positive training set of high-quality' ' variants for VQSR'),
    },
    'NEGATIVE_TRAIN_SITE': {
        'Description': ('Variant was used to build the negative training set of low-quality' ' variants for VQSR'),
    },
    'positive_train_site': {
        'Description': ('Variant was used to build the positive training set of high-quality' ' variants for VQSR'),
    },
    'negative_train_site': {
        'Description': ('Variant was used to build the negative training set of low-quality' ' variants for VQSR'),
    },
    'BaseQRankSum': {
        'Description': ('Z-score from Wilcoxon rank sum test of alternate vs. reference base' ' qualities'),
    },
    'VarDP': {
        'Description': ('Depth over variant genotypes (does not include depth of reference samples)'),
    },
    'VQSLOD': {
        'Description': (
            'Log-odds ratio of being a true variant versus being a false positive under'
            ' the trained VQSR Gaussian mixture model'
        ),
    },
    'culprit': {
        'Description': 'Worst-performing annotation in the VQSR Gaussian mixture model',
    },
    'decoy': {'Description': 'Variant falls within a reference decoy region'},
    'lcr': {'Description': 'Variant falls within a low complexity region'},
    'non_par': {
        'Description': ('Variant (on sex chromosome) falls outside a pseudoautosomal region'),
    },
    'segdup': {'Description': 'Variant falls within a segmental duplication region'},
    'outside_capture_region': {
        'Description': 'Variant falls outside of the OurDNA exome capture regions.',
    },
    'outside_calling_region': {
        'Description': ('Variant falls outside of the OurDNA exome capture regions plus 150 bp padding.'),
    },
    'rf_positive_label': {
        'Description': ('Variant was labelled as a positive example for training of random forest' ' model'),
    },
    'rf_negative_label': {
        'Description': ('Variant was labelled as a negative example for training of random forest' ' model'),
    },
    'rf_label': {'Description': 'Random forest training label'},
    'rf_train': {'Description': 'Variant was used in training random forest model'},
    'rf_tp_probability': {
        'Description': ('Probability of a called variant being a true variant as determined by' ' random forest model'),
    },
    'transmitted_singleton': {
        'Description': (
            'Variant was a callset-wide doubleton that was transmitted within a family'
            ' from a parent to a child (i.e., a singleton amongst unrelated samples in'
            ' cohort)'
        ),
    },
    'sibling_singleton': {
        'Description': (
            'Variant was a callset-wide doubleton that was present only in two siblings'
            ' (i.e., a singleton amongst unrelated samples in cohort).'
        ),
    },
    'original_alleles': {'Description': 'Alleles before splitting multiallelics'},
    'variant_type': {
        'Description': 'Variant type (snv, indel, multi-snv, multi-indel, or mixed)',
    },
    'allele_type': {
        'Description': 'Allele type (snv, insertion, deletion, or mixed)',
    },
    'n_alt_alleles': {
        'Number': '1',
        'Description': 'Total number of alternate alleles observed at variant locus',
    },
    'was_mixed': {'Description': 'Variant type was mixed'},
    'has_star': {
        'Description': (
            'Variant locus coincides with a spanning deletion (represented by a star)'
            ' observed elsewhere in the callset'
        ),
    },
    'AS_pab_max': {
        'Number': 'A',
        'Description': (
            'Maximum p-value over callset for binomial test of observed allele balance'
            ' for a heterozygous genotype, given expectation of 0.5'
        ),
    },
    'monoallelic': {
        'Description': 'All samples are homozygous alternate for the variant',
    },
    'only_het': {'Description': 'All samples are heterozygous for the variant'},
    'QUALapprox': {
        'Number': '1',
        'Description': 'Sum of PL[0] values; used to approximate the QUAL score',
    },
    'AS_SB_TABLE': {
        'Number': '.',
        'Description': ('Allele-specific forward/reverse read counts for strand bias tests'),
    },
    # Hard coding this one due to `grpmax` not being a genetic ancestry to iterate over
    'faf95_grpmax': {
        'Number': 'A',
        'Description': 'Filtering allele frequency (using Poisson 95% CI) for the genetic ancestry group with the maximum allele',
    },
}


def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.Table:
    """
    Create a Table ready for vcf export.

    In particular, the following conversions are done:
        - All int64 are coerced to int32
        - Fields specified by `pipe_delimited_annotations` are converted from arrays to pipe-delimited strings

    :param ht: Input Table.
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct).
    :return: Table ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), '')), '|')

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by('locus', 'alleles')

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype('int64'):
            logger.warning(
                f'Coercing field info.{f} from int64 to int32 for VCF output. Value'
                ' will be capped at int32 max value.',
            )
            info_type_convert_expr.update(
                {
                    f: hl.or_missing(
                        hl.is_defined(ht.info[f]),
                        hl.int32(hl.min(2**31 - 1, ht.info[f])),
                    ),
                },
            )
        elif ft == hl.dtype('array<int64>'):
            logger.warning(
                f'Coercing field info.{f} from array<int64> to array<int32> for VCF'
                f' output. Array values will be capped at int32 max value.',
            )
            info_type_convert_expr.update(
                {
                    f: ht.info[f].map(
                        lambda x: hl.or_missing(
                            hl.is_defined(x),
                            hl.int32(hl.min(2**31 - 1, x)),
                        ),
                    ),
                },
            )

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != 'AS_SB_TABLE':
            info_expr[f] = '|' + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if 'SB' in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr['SB'] = ht.info.SB[0].extend(ht.info.SB[1])

    if 'AS_SB_TABLE' in ht.info:
        info_expr['AS_SB_TABLE'] = get_pipe_expr(
            ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ',')),
        )

    # Annotate with new expression
    ht = ht.annotate(info=ht.info.annotate(**info_expr))

    return ht


def format_validated_ht_for_export(
    ht: hl.Table,
    data_type: str = 'exomes',
    vcf_info_reorder: List[str] = VCF_INFO_REORDER,
    info_fields_to_drop: Optional[List[str]] = None,
) -> Tuple[hl.Table, List[str]]:
    """
    Format validated HT for export.

    Drop downsamplings frequency stats from info, rearrange info, and make sure fields
    are VCF compatible.

    :param ht: Validated HT.
    :param data_type: Data type to format validated HT for. One of "exomes" or "genomes".
        Default is "exomes".
    :param vcf_info_reorder: Order of VCF INFO fields. These will be placed in front of
        all other fields in the order specified.
    :param info_fields_to_drop: List of info fields to drop from the info struct.
    :return: Formatted HT and list rename row annotations.
    """
    if info_fields_to_drop is None:
        info_fields_to_drop = []
    logger.info('Dropping fields from info struct...')

    if data_type == 'joint':
        for dt in ['exomes', 'genomes', 'joint']:
            info_fields_to_drop.extend(
                [
                    f'age_hist_het_bin_edges_{dt}',
                    f'age_hist_hom_bin_edges_{dt}',
                ],
            )
    else:
        logger.info('Add age_histogram bin edges to info fields to drop...')
        info_fields_to_drop.extend(['age_hist_het_bin_edges', 'age_hist_hom_bin_edges'])

        logger.info('Adding "SB" to info fields to drop...')
        info_fields_to_drop.append('SB')

    logger.info('Dropping the following fields from info struct:')
    pprint(info_fields_to_drop)
    ht = ht.annotate(info=ht.info.drop(*info_fields_to_drop))

    logger.info('Dropping _"adj" from info fields...')
    row_annots = list(ht.info)
    new_row_annots = [x.replace('_adj', '') for x in row_annots]
    info_annot_mapping = dict(
        zip(new_row_annots, [ht.info[f'{x}'] for x in row_annots]),
    )
    ht = ht.transmute(info=hl.struct(**info_annot_mapping))

    logger.info('Adjusting VCF incompatible types...')
    if data_type != 'joint':
        # Reformat AS_SB_TABLE for use in adjust_vcf_incompatible_types
        ht = ht.annotate(
            info=ht.info.annotate(
                AS_SB_TABLE=hl.array([ht.info.AS_SB_TABLE[:2], ht.info.AS_SB_TABLE[2:]]),
            ),
        )
    # The Table is already split so there are no annotations that need to be
    # pipe delimited.
    ht = adjust_vcf_incompatible_types(ht, pipe_delimited_annotations=[])

    logger.info('Rearranging fields to desired order...')
    if data_type == 'joint':
        special_items = {'exomes': 'exomes_filters', 'genomes': 'genomes_filters'}
        new_vcf_info_reorder = []
        for dt in ['joint', 'exomes', 'genomes']:
            new_vcf_info_reorder += [f'{f}_{dt}' for f in vcf_info_reorder]
            if dt in special_items:
                new_vcf_info_reorder.append(special_items[dt])
        vcf_info_reorder = new_vcf_info_reorder

    ht = ht.annotate(
        info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder)),
    )
    return ht, new_row_annots


def select_type_from_joint_ht(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Select all fields from the joint HT that are relevant to `data_type`.

    :param ht: Joint release HT.
    :param data_type: Data type to select in joint HT. One of "exomes", "genomes", or
        "joint".
    :return: Joint HT with fields relevant to `data_type`.
    """
    global_fields = [f'{data_type}_globals']
    row_fields = [data_type, 'region_flags']
    if data_type == 'joint':
        row_fields.append('freq_comparison_stats')
    ht = ht.select_globals(*global_fields)
    ht = ht.select(*row_fields)
    ht = ht.transmute_globals(**ht[f'{data_type}_globals'])
    ht = ht.transmute(**ht[data_type])
    return ht


def unfurl_nested_annotations(
    ht: hl.Table,
    entries_to_remove: Set[str] | None = None,
    data_type: str = 'exomes',
    joint_included: bool = False,
    hist_prefix: str = '',
    freq_comparison_included: bool = False,
    for_joint_validation: bool = False,
) -> tuple[hl.expr.StructExpression, Set[str] | None, Dict[str, str]]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays.

    The values of the returned dictionary are Hail Expressions describing how to access
    the corresponding values.

    :param ht: Table containing the nested variant annotation arrays to be unfurled.
    :param entries_to_remove: Optional Set of frequency entries to remove for vcf_export.
    :param data_type: Data type to unfurl nested annotations for. One of "exomes",
        "genomes", or "joint". Default is "exomes".
    :param joint_included: Whether joint frequency data is included in the exome or
        genome HT. Default is False.
    :param hist_prefix: Prefix to use for histograms. Default is "".
    :param freq_comparison_included: Whether frequency comparison data is included in
        the HT. Default is False.
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is
        False.
    :return: StructExpression containing variant annotations and their corresponding
        expressions and updated entries, set of frequency entries to remove from the
        VCF, and a dict of fields to rename when `for_joint_validation` is True.
    """
    expr_dict = {}
    rename_dict = {}

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=afr_XX, i=31)
    logger.info('Unfurling freq data...')
    freq_idx = hl.eval(ht.freq_index_dict)
    for k, i in freq_idx.items():
        for f in ht.freq[0].keys():
            field_name = f if f != 'homozygote_count' else 'nhomalt'
            expr_dict[f'{field_name}_{k}'] = ht.freq[i][f]
            if for_joint_validation:
                rename_dict[f'{field_name}_{k}'] = f'{field_name}_{data_type}_{k}'

    if joint_included:
        logger.info('Unfurling joint freq data...')
        joint_freq_idx = hl.eval(ht.joint_freq_index_dict)
        expr_dict.update(
            {
                f'{f if f != "homozygote_count" else "nhomalt"}_joint_{k}': (ht.joint_freq[i][f])
                for k, i in joint_freq_idx.items()
                for f in ht.joint_freq[0].keys()
            },
        )

    # This creates fields like grpmax, AC_grpmax_non_ukb...
    logger.info('Adding grpmax data...')
    # Our tables use popmax (exomes/genomes) vs grpmax (joint) with 'pop' vs 'gen_anc' keys respectively.
    # Our pipeline uses pop_max_expr() with {'group', 'pop'} keys, while gnomAD v4 uses the equivalent
    # grpmax_expr() with {'group', 'gen_anc'} keys.
    grpmax_idx = ht.grpmax if for_joint_validation else ht.popmax
    grpmax_dict = {'grpmax': grpmax_idx.pop}
    grpmax_rename = {f: f if f != 'homozygote_count' else 'nhomalt' for f in grpmax_idx.keys() if f != 'pop'}
    grpmax_dict.update(
        {f'{v}_grpmax': grpmax_idx[k] for k, v in grpmax_rename.items()},
    )
    if for_joint_validation:
        rename_dict['grpmax'] = f'grpmax_{data_type}'
        rename_dict.update(
            {f'{v}_grpmax': f'{v}_grpmax_{data_type}' for v in grpmax_rename.values()},
        )
    expr_dict.update(grpmax_dict)

    if joint_included:
        logger.info('Adding joint grpmax data...')
        joint_grpmax_idx = ht.joint_grpmax
        joint_grpmax_dict = {'grpmax_joint': joint_grpmax_idx.pop}
        joint_grpmax_dict.update(
            {
                f'{f if f != "homozygote_count" else "nhomalt"}_grpmax_joint': (joint_grpmax_idx[f])
                for f in [f for f in joint_grpmax_idx._fields if f != 'pop']
            },
        )
        expr_dict.update(joint_grpmax_dict)

    logger.info('Unfurling faf data...')
    faf_idx = hl.eval(ht.faf_index_dict)
    for k, i in faf_idx.items():
        for f in ht.faf[0].keys():
            expr_dict[f"{f}_{k}"] = ht.faf[i][f]
            if for_joint_validation:
                rename_dict[f"{f}_{k}"] = f"{f}_{data_type}_{k}"

    logger.info('Unfurling fafmax data...')
    fafmax_idx = ht.fafmax

    fafmax_dict = {f'fafmax_{f}': fafmax_idx[f] for f in fafmax_idx.keys()}
    if for_joint_validation:
        rename_dict.update(
            {f'fafmax_{f}': f'fafmax_{f}_{data_type}' for f in fafmax_idx.keys()},
        )
    expr_dict.update(fafmax_dict)

    if joint_included:
        logger.info('Unfurling joint faf data...')
        joint_faf_idx = hl.eval(ht.joint_faf_index_dict)
        expr_dict.update(
            {f'{f}_joint_{k}': ht.joint_faf[i][f] for f in ht.joint_faf[0].keys() for k, i in joint_faf_idx.items()},
        )

        logger.info('Unfurling joint fafmax data...')
        joint_fafmax_idx = ht.joint_fafmax
        joint_fafmax_dict = {
            f'fafmax_{f if f != "joint_fafmax_data_type" else "data_type"}_joint': (joint_fafmax_idx[f])
            for f in joint_fafmax_idx.keys()
        }
        expr_dict.update(joint_fafmax_dict)

    logger.info('Unfurling age hists...')
    age_hists = ['age_hist_het', 'age_hist_hom']
    hist_idx = ht.histograms.age_hists if for_joint_validation else ht
    for hist in age_hists:
        for f in hist_idx[hist].keys():
            expr_dict[f'{hist}_{f}'] = hl.delimit(hist_idx[hist][f], delimiter='|') if 'bin' in f else hist_idx[hist][f]
            if for_joint_validation:
                rename_dict[f'{hist}_{f}'] = f'{hist}_{f}_{data_type}'

    logger.info('Unfurling variant quality histograms...')
    # Add underscore to hist_prefix if it isn't empty
    if hist_prefix != '':
        hist_prefix += '_'

    # Histograms to export are:
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    # We previously dropped:
    # _n_smaller for all hists
    # _bin_edges for all hists
    # _n_larger for all hists EXCEPT DP hists
    for hist in HISTS:
        hist_dict = {
            f'{hist}_bin_freq': hl.delimit(
                ht.histograms.qual_hists[hist].bin_freq,
                delimiter='|',
            ),
        }
        expr_dict.update(hist_dict)
        if for_joint_validation:
            rename_dict.update(
                {f'{hist}_bin_freq': f'{hist_prefix}{hist}_bin_freq_{data_type}'},
            )

        if 'dp' in hist:
            expr_dict.update(
                {f'{hist}_n_larger': ht.histograms.qual_hists[hist].n_larger},
            )
            if for_joint_validation:
                rename_dict.update(
                    {f'{hist}_n_larger': f'{hist_prefix}{hist}_n_larger_{data_type}'},
                )

    if freq_comparison_included:
        logger.info('Unfurling contingency table test results...')
        contingency_idx = hl.eval(ht.freq_index_dict)
        for k, i in contingency_idx.items():
            for f in ht.freq_comparison_stats.contingency_table_test[0].keys():
                key = f'CTT_{f}_{k}'
                expr = ht.freq_comparison_stats.contingency_table_test[i][f]
                expr_dict[key] = expr

        logger.info('Unfurling Cochran-Mantel-Haenszel test results...')
        expr_dict['CMH_chisq'] = ht.freq_comparison_stats.cochran_mantel_haenszel_test.chisq
        expr_dict['CMH_p_value'] = ht.freq_comparison_stats.cochran_mantel_haenszel_test.p_value
        logger.info('Unfurling unionized stats...')
        expr_dict['stat_union_p_value'] = ht.freq_comparison_stats.stat_union.p_value
        expr_dict['stat_union_test_name'] = ht.freq_comparison_stats.stat_union.stat_test_name
        expr_dict['stat_union_gen_ancs'] = ht.freq_comparison_stats.stat_union.gen_ancs

    return hl.struct(**expr_dict), entries_to_remove, rename_dict


def make_info_expr(
    t: hl.Table,
    data_type: str = 'exomes',
    for_joint_validation: bool = False,
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param t: Table containing variant annotations to be reformatted for VCF export.
    :param data_type: Data type to make info expression for. One of "exomes", "genomes",
        or "joint". Default is "exomes".
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is False.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    """
    vcf_info_dict = {}

    # Set data type to joint if for_joint_validation is True so the correct region flag
    # fields are used.
    if for_joint_validation:
        data_type = 'joint'

    if 'region_flags' in t.row:
        # Add region_flag to info dict
        for field in REGION_FLAG_FIELDS[data_type]:
            vcf_info_dict[field] = t['region_flags'][f'{field}']

    if for_joint_validation:
        return vcf_info_dict

    # Add site-level annotations and AS annotations to vcf_info_dict
    for field in SITE_FIELDS[data_type] + AS_FIELDS:
        vcf_info_dict[field] = t['release_ht_info'][f'{field}']

    for field in AS_VQSR_FIELDS:
        # NOTE: VQSR results are nested in the info struct (and not 'vqsr_results') and are also
        # quasi-AS not true-AS in our release tables.
        # We also moved the info fields to `release_ht_info` field.
        vcf_info_dict[field] = t.release_ht_info[f'{field}']

    # Add allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS[data_type]:
        vcf_info_dict[field] = t['allele_info'][f'{field}']

    # Add vep annotations to info dict
    vcf_info_dict['vep'] = t['vep']

    # Add monoallelic field to info dict
    vcf_info_dict['monoallelic'] = t['monoallelic']

    return vcf_info_dict


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression,
    csq_fields: str = VEP_CSQ_FIELDS[CURRENT_VEP_VERSION],
    has_polyphen_sift: bool = True,
) -> hl.expr.ArrayExpression:
    """
    Given a VEP Struct, returns and array of VEP VCF CSQ strings (one per consequence in the struct).

    The fields and their order will correspond to those passed in `csq_fields`, which corresponds to the
    VCF header that is required to interpret the VCF CSQ INFO field.

    Note that the order is flexible and that all fields that are in the default value are supported.
    These fields will be formatted in the same way that their VEP CSQ counterparts are.

    While other fields can be added if their name are the same as those in the struct. Their value will be the result of calling
    hl.str(), so it may differ from their usual VEP CSQ representation.

    :param vep_expr: The input VEP Struct
    :param csq_fields: The | delimited list of fields to include in the CSQ (in that order), default is the CSQ fields of the CURRENT_VEP_VERSION.
    :param has_polyphen_sift: Whether the input VEP Struct has PolyPhen and SIFT annotations. Default is True.
    :return: The corresponding CSQ strings
    """
    _csq_fields = [f.lower() for f in csq_fields.split('|')]

    def get_csq_from_struct(
        element: hl.expr.StructExpression,
        feature_type: str,
    ) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                'allele': element.variant_allele,
                'consequence': hl.delimit(element.consequence_terms, delimiter='&'),
                'feature_type': feature_type,
                'feature': (
                    element.transcript_id
                    if 'transcript_id' in element
                    else (
                        element.regulatory_feature_id
                        if 'regulatory_feature_id' in element
                        else (element.motif_feature_id if 'motif_feature_id' in element else '')
                    )
                ),
                'variant_class': vep_expr.variant_class,
            },
        )

        # Add exception for transcripts
        if feature_type == 'Transcript':
            transcript_dict = {
                'canonical': hl.if_else(element.canonical == 1, 'YES', ''),
                'ensp': element.protein_id,
                'gene': element.gene_id,
                'symbol': element.gene_symbol,
                'symbol_source': element.gene_symbol_source,
                'cdna_position': hl.str(element.cdna_start)
                + hl.if_else(
                    element.cdna_start == element.cdna_end,
                    '',
                    '-' + hl.str(element.cdna_end),
                ),
                'cds_position': hl.str(element.cds_start)
                + hl.if_else(
                    element.cds_start == element.cds_end,
                    '',
                    '-' + hl.str(element.cds_end),
                ),
                'mirna': hl.delimit(element.mirna, '&') if 'mirna' in element else None,
                'protein_position': hl.str(element.protein_start)
                + hl.if_else(
                    element.protein_start == element.protein_end,
                    '',
                    '-' + hl.str(element.protein_end),
                ),
                'uniprot_isoform': (hl.delimit(element.uniprot_isoform, '&') if 'uniprot_isoform' in element else None),
            }
            # Retain transcript dict updates only for fields that exist in the csq
            # fields.
            transcript_dict = {
                k: v for k, v in transcript_dict.items() if k in [x.lower() for x in csq_fields.split('|')]
            }
            fields.update(transcript_dict)

            if has_polyphen_sift:
                fields.update(
                    {
                        'sift': (element.sift_prediction + '(' + hl.format('%.3f', element.sift_score) + ')'),
                        'polyphen': (
                            element.polyphen_prediction + '(' + hl.format('%.3f', element.polyphen_score) + ')'
                        ),
                    },
                )
            fields.update(
                {
                    'domains': hl.delimit(
                        element.domains.map(lambda d: d.db + ':' + d.name),
                        '&',
                    ),
                },
            )
        elif feature_type == 'MotifFeature':
            fields['motif_score_change'] = hl.format('%.3f', element.motif_score_change)
            if 'transcription_factors' in element:
                fields['transcription_factors'] = hl.delimit(
                    element.transcription_factors,
                    '&',
                )

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, '')), '') for f in _csq_fields],
            '|',
        )

    csq = hl.empty_array(hl.tstr)
    for feature_field, feature_type in [
        ('transcript_consequences', 'Transcript'),
        ('regulatory_feature_consequences', 'RegulatoryFeature'),
        ('motif_feature_consequences', 'MotifFeature'),
        ('intergenic_consequences', 'Intergenic'),
    ]:
        csq = csq.extend(
            hl.or_else(
                vep_expr[feature_field].map(
                    lambda x: get_csq_from_struct(x, feature_type=feature_type),
                ),
                hl.empty_array(hl.tstr),
            ),
        )

    return hl.or_missing(hl.len(csq) > 0, csq)


def process_vep_csq_header(vep_csq_header: str = VEP_CSQ_HEADER) -> str:
    """
    Process VEP CSQ header string, delimited by '|', to remove polyphen and sift annotations.

    :param vep_csq_header: VEP CSQ header.
    :return: Processed VEP CSQ header.
    """
    logger.info('Processing VEP CSQ header...')
    csq_fields = vep_csq_header.split('|')
    csq_fields = [f for f in csq_fields if f not in ['PolyPhen', 'SIFT']]
    vep_csq_header = '|'.join(csq_fields)
    return vep_csq_header


def _freq(ds, *args, **kwargs):
    return ds.freq[ds.freq_index_dict[_freq_index_key(*args, **kwargs)]]


def get_filters_expr(ht: hl.Table, score_cutoffs: dict) -> hl.Table:
    inbreeding_coeff_cutoff = config_retrieve(['large_cohort', 'browser', 'inbreeding_coeff_cutoff'])

    ac = _freq(ht, subset=None).AC
    filters = {
        'inbreeding_coeff': ht.inbreeding_coeff[0] < inbreeding_coeff_cutoff,
        'AC0': ac == 0,
        'AS_lowqual': ht.AS_lowqual,
        'AS_VQSR': hl.is_missing(ht.info['AS_VQSLOD']),
    }
    # NOTE: `score_cutoffs` are accessed during browser prep, need to read in stage output.
    snv_indel_expr = {'snv': hl.is_snp(ht.alleles[0], ht.alleles[1])}
    snv_indel_expr['indel'] = ~snv_indel_expr['snv']
    if score_cutoffs is not None:
        for var_type, score_cut in score_cutoffs.items():
            filters['AS_VQSR'] = filters['AS_VQSR'] | (
                snv_indel_expr[var_type] & (ht.info.AS_VQSLOD < score_cut.min_score)
            )
    return ht.annotate(filters=add_filters_expr(filters=filters))


def prepare_ht_for_validation(
    ht: hl.Table,
    data_type: str = 'exomes',
    freq_entries_to_remove: Optional[Set[str]] = None,
    vcf_info_reorder: Optional[List[str]] = None,
    joint_included: bool = False,
    for_joint_validation: bool = True,
    freq_comparison_included: bool = False,
    score_cutoffs: dict[str, hl.expr.StructExpression] | None = None,
) -> hl.Table:
    """
    Prepare HT for validity checks and export.

    :param ht: Release Hail Table.
    :param data_type: Data type to prepare HT for. One of "exomes", "genomes", or
        "joint". Default is "exomes".
    :param freq_entries_to_remove: List of entries to remove from freq.
    :param vcf_info_reorder: Order of VCF INFO fields.
    :param joint_included: Whether joint frequency data is included in the HT. Default
        is False.
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is
        True.
    :param freq_comparison_included: Whether frequency comparison data is included in
        the HT. Default is False.
    :param score_cutoffs: Dictionary of score cutoffs for VQSR filtering. Keys are
        "snv" and "indel" and values are StructExpressions with fields "min_score" and
        "bin_id".
    :return: Hail Table prepared for validity checks and export and a dictionary of
        fields to rename when `for_joint_validation` is True.
    """
    logger.info(
        'Unfurling nested gnomAD frequency annotations and add to INFO field...',
    )
    info_struct, freq_entries_to_remove, rename_dict = unfurl_nested_annotations(
        ht,
        entries_to_remove=freq_entries_to_remove,
        data_type=data_type,
        joint_included=joint_included,
        freq_comparison_included=freq_comparison_included,
        for_joint_validation=for_joint_validation,
    )

    logger.info('Constructing INFO field')
    # Remove SIFT and Polyphen from CSQ fields or they will be inserted with
    # missing values by vep_struct_to_csq. These fields are processed separately
    # as in silico annotations.
    csq_fields = '|'.join(
        [c for c in VEP_CSQ_FIELDS[CURRENT_VEP_VERSION].split('|') if c != 'PolyPhen' and c != 'SIFT'],
    )

    if for_joint_validation:
        ann_expr = {'info': info_struct}
        if 'region_flags' in ht.row:
            ann_expr['region_flags'] = ht.region_flags
    else:
        ann_expr = {
            'region_flag': ht.region_flags,
            'release_ht_info': ht.info,
            'info': info_struct,
            'rsid': ht.rsid,  # NOTE: In exomes/genomes tables rsid field is not an array
            'vep': vep_struct_to_csq(
                ht.vep,
                csq_fields=csq_fields,
                has_polyphen_sift=False,
            ),
        }

    ht = ht.annotate(**ann_expr)

    # Add variant annotations to INFO field
    # This adds the following:
    #   region flag for problematic regions
    #   annotations in ht.release_ht_info (site and allele-specific annotations),
    #   info struct (unfurled data obtained above),
    #   dbSNP rsIDs
    #   all VEP annotations
    ht = ht.annotate(
        info=ht.info.annotate(
            **make_info_expr(
                ht,
                data_type=data_type,
                for_joint_validation=for_joint_validation,
            ),
        ),
    )

    if for_joint_validation:
        ht = ht.annotate_globals(
            freq_entries_to_remove=(freq_entries_to_remove if freq_entries_to_remove else hl.empty_set(hl.tstr)),
        )
    else:
        ht = ht.annotate_globals(
            vep_csq_header=process_vep_csq_header(VEP_CSQ_HEADER),
            freq_entries_to_remove=(freq_entries_to_remove if freq_entries_to_remove else hl.empty_set(hl.tstr)),
        )

    # Select relevant fields for VCF export
    if for_joint_validation:
        if 'filters' in ht.row:
            filters_expr = ht.filters
        else:
            filters_expr = hl.empty_set(hl.tstr)
        ht = ht.select('info', filters=filters_expr)
    else:
        # Our frequencies tables do not have a `filters` field. This is only annotated during browser prep. Doing it here.
        inbreeding_coeff_cutoff = config_retrieve(['large_cohort', 'browser', 'inbreeding_coeff_cutoff'])
        ac = _freq(ht, subset=None).AC
        filters = {
            'inbreeding_coeff': ht.inbreeding_coeff[0] < inbreeding_coeff_cutoff,
            'AC0': ac == 0,
            'AS_lowqual': ht.AS_lowqual,
            'AS_VQSR': hl.is_missing(ht.info['AS_VQSLOD']),
        }
        # NOTE: `score_cutoffs` are accessed during browser prep, need to read in stage output.
        snv_indel_expr = {'snv': hl.is_snp(ht.alleles[0], ht.alleles[1])}
        snv_indel_expr['indel'] = ~snv_indel_expr['snv']
        if score_cutoffs is not None:
            for var_type, score_cut in score_cutoffs.items():
                filters['AS_VQSR'] = filters['AS_VQSR'] | (
                    snv_indel_expr[var_type] & (ht.info.AS_VQSLOD < score_cut.min_score)
                )
        ht = ht.annotate(filters=add_filters_expr(filters=filters))

        ht = ht.select('info', 'filters', 'rsid')

    if vcf_info_reorder:
        logger.info('Rearranging fields to desired order...')
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder)),
        )

    return ht, rename_dict


def get_joint_filters(ht: hl.Table) -> hl.Table:
    """
    Transform exomes and genomes filters to joint filters.

    :param ht: Input Table.
    :return: Table with joint filters transformed from exomes and genomes filters.
    """
    # NOTE: I'm not sure if gnomAD uses AC0 as a filter in the joint release. If AC0 is
    # used in genomes or exomes, then it should be added to the logic below. Currently,
    # AC0 is not included in the joint release filters and should possibly be a separate
    # item in the filters list based on ac = ht.joint.grpmax.AC
    # Or is it ht.joint.freq.AC?
    exomes_filters = ht.info.exomes_filters
    genomes_filters = ht.info.genomes_filters
    ht = ht.annotate(
        filters=hl.set(
            hl.case()
            .when(
                ((hl.len(exomes_filters) == 0) | hl.is_missing(exomes_filters))
                & ((hl.len(genomes_filters) == 0) | hl.is_missing(genomes_filters)),
                ['PASS'],
            )
            .when(
                (hl.len(exomes_filters) != 0) & ((hl.len(genomes_filters) == 0) | hl.is_missing(genomes_filters)),
                ['EXOMES_FILTERED'],
            )
            .when(
                ((hl.len(exomes_filters) == 0) | hl.is_missing(exomes_filters)) & (hl.len(genomes_filters) != 0),
                ['GENOMES_FILTERED'],
            )
            .when(
                (hl.len(exomes_filters) != 0) & (hl.len(genomes_filters) != 0),
                ['BOTH_FILTERED'],
            )
            .default(['MISSING_FILTERS']),
        ),
    )
    return ht


def make_vcf_filter_dict(
    snp_cutoff: Optional[float] = None,
    indel_cutoff: Optional[float] = None,
    inbreeding_cutoff: Optional[float] = None,
    variant_qc_filter: str = 'RF',
    joint: bool = False,
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    Generates descriptions for:
        - AC0 filter
        - InbreedingCoeff filter
        - Variant QC filter (RF or AS_VQSR)
        - PASS (passed all variant filters)

    :param snp_cutoff: Minimum SNP cutoff score from random forest model.
    :param indel_cutoff: Minimum indel cutoff score from random forest model.
    :param inbreeding_cutoff: Inbreeding coefficient hard cutoff.
    :param variant_qc_filter: Method used for variant QC filter. One of 'RF' or 'AS_VQSR'. Default is 'RF'.
    :param joint: Whether the filter dictionary is for the joint release. Default is False.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
    """
    variant_qc_filter_dict = {
        'RF': {
            'Description': (
                f'Failed random forest filtering thresholds of {snp_cutoff} for SNPs'
                f' and {indel_cutoff} for indels (probabilities of being a true'
                ' positive variant)'
            ),
        },
        'AS_VQSR': {
            'Description': (
                f'Failed VQSR filtering thresholds of {snp_cutoff} for SNPs and' f' {indel_cutoff} for indels'
            ),
        },
    }

    if variant_qc_filter not in variant_qc_filter_dict:
        raise ValueError(
            f'{variant_qc_filter} is not a valid value for "variant_qc_filter". It must be "RF" or "AS_VQSR"',
        )
    if not joint and (snp_cutoff is None or indel_cutoff is None or inbreeding_cutoff is None):
        raise ValueError(
            'snp_cutoff, indel_cutoff, and inbreeding_cutoff must be specified to generate filter descriptions.',
        )

    if joint:
        filter_dict = {
            'PASS': {
                'Description': 'Either passed all variant filters in both exomes and '
                'genomes, or passed all variant filters in either '
                'exomes or genomes while being absent from the other '
                'dataset',
            },
            'EXOMES_FILTERED': {
                'Description': 'Failed variant filters in the exomes dataset and either '
                'passed all variant filters in the genomes dataset or the variant was '
                'not present in the genomes dataset. Refer to "exomes_filters" within '
                'INFO for more information',
            },
            'GENOMES_FILTERED': {
                'Description': 'Failed variant filters in the genomes dataset and either '
                'passed all variant filters in the exomes dataset or the variant was '
                'not present in the exomes dataset. Refer to "genomes_filters" within '
                'INFO for more information',
            },
            'BOTH_FILTERED': {
                'Description': 'Failed variant filters in both exomes and genomes datasets. '
                'Refer to "exomes_filters" and "genomes_filters" within INFO for more information',
            },
        }
    else:
        # NOTE: Are these accurate?
        filter_dict = {
            'AC0': {
                'Description': (
                    'Allele count is zero after filtering out low-confidence genotypes (GQ'
                    ' < 20; DP < 10; and AB < 0.2 for het calls)'
                ),
            },
            'inbreeding_coeff': {
                'Description': f'Inbreeding coefficient < {inbreeding_cutoff}',
            },
            'AS_lowqual': {
                'Description': (
                    'Whether the variant falls below a low quality threshold and was '
                'excluded from the OurDNA dataset'
                ),
            },
            'PASS': {'Description': 'Passed all variant filters'},
            variant_qc_filter: variant_qc_filter_dict[variant_qc_filter],
        }

    return filter_dict


def prepare_vcf_header_dict(
    ht: hl.Table,
    validated_ht: Optional[hl.Table],
    info_fields: List[str],
    bin_edges: Dict[str, str],
    age_hist_distribution: str,
    subset_list: List[str],
    pops: Dict[str, str],
    data_type: str = 'exomes',
    joint_included: bool = False,
    freq_comparison_included: bool = False,
    extra_suffix: str | None = None,
    extra_description_text: str | None = None,
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param ht: Input Table
    :param validated_ht: Validated HT with unfurled info fields.
    :param info_fields: List of info fields to add to the info dict.
    :param bin_edges: Dictionary of variant annotation histograms and their associated
        bin edges.
    :param age_hist_distribution: Pipe-delimited string of overall age histogram bin
        frequency.
    :param subset_list: List of sample subsets in dataset.
    :param pops: List of sample global genetic ancestry group names for gnomAD data type.
    :param data_type: Data type to prepare VCF header for. One of "exomes" or "genomes".
        Default is "exomes".
    :param joint_included: Whether joint frequency data is included in the HT. Default is False.
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :param extra_suffix: Suffix to add to INFO field.
    :param extra_description_text: Extra description text to add to INFO field.
    :return: Prepared VCF header dictionary.
    """
    if data_type != 'joint':
        logger.info('Making FILTER dict for VCF...')
        filter_dict = make_vcf_filter_dict(
            hl.eval(ht.filtering_model.snv_cutoff.min_score),
            hl.eval(ht.filtering_model.indel_cutoff.min_score),
            inbreeding_cutoff=hl.eval(ht.inbreeding_coeff_cutoff),
            variant_qc_filter=hl.eval(ht.filtering_model.filter_name),
        )
        # subset = '' represents full dataset in VCF header construction, the
        # logic in gnomad_methods is built around this.
        subset_list.extend(['', 'joint'] if joint_included else [''])

    logger.info('Making INFO dict for VCF...')
    vcf_info_dict = populate_info_dict(
        info_fields=info_fields,
        bin_edges=bin_edges,
        age_hist_distribution=age_hist_distribution,
        subset_list=subset_list,
        pops=pops,
        data_type=data_type,
        freq_comparison_included=freq_comparison_included,
        extra_suffix=extra_suffix,
        extra_description_text=extra_description_text,
    )

    if data_type != 'joint':
        if validated_ht is not None:
            vcf_info_dict.update(
                # NOTE: Check that the consequences in validated_ht.vep_csq_header are all in our dataset!
                {'vep': {'Description': hl.eval(validated_ht.vep_csq_header)}},
            )

    # Adjust keys to remove adj tags before exporting to VCF.
    new_vcf_info_dict = {i.replace('_adj', ''): j for i, j in vcf_info_dict.items()}

    if data_type == 'joint':
        header_dict = new_vcf_info_dict
    else:
        header_dict = {
            'info': new_vcf_info_dict,  # type: ignore[dict-item]
            'filter': filter_dict,  # type: ignore[dict-item]
        }

    return header_dict


def create_label_groups(
    gen_ancs: List[str],
    sexes: List[str] = SEXES,
    all_groups: List[str] = GROUPS,
    gen_anc_sex_groups: List[str] = ['adj'],
) -> List[Dict[str, List[str]]]:
    """
    Generate a list of label group dictionaries needed to populate info dictionary.

    Label dictionaries are passed as input to `make_info_dict`.

    :param gen_ancs: List of genetic ancestry group names.
    :param sexes: List of sample sexes.
    :param all_groups: List of data types (raw, adj). Default is `GROUPS`, which is ["raw", "adj"].
    :param gen_anc_sex_groups: List of data types (raw, adj) to populate with gen_ancs and sexes. Default is ["adj"].
    :return: List of label group dictionaries.
    """
    return [
        # This is to capture raw frequency fields, which are
        # not stratified by sex or genetic ancestry group (e.g., only AC_raw
        # exists, not AC_XX_raw)
        dict(group=all_groups),
        dict(group=gen_anc_sex_groups, sex=sexes),
        dict(group=gen_anc_sex_groups, gen_anc=gen_ancs),
        dict(group=gen_anc_sex_groups, gen_anc=gen_ancs, sex=sexes),
    ]


def populate_subset_info_dict(
    subset: str,
    description_text: str,
    data_type: str = 'exomes',
    pops: Dict[str, str] = GEN_ANC_GROUPS['exomes'],  # type: ignore[assignment]
    faf_pops: Dict[str, List[str]] = FAF_GEN_ANC_GROUPS,
    sexes: List[str] = SEXES,
    label_delimiter: str = '_',
    freq_comparison_included: bool = False,
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` to populate INFO dictionary for the requested `subset`.

    Creates:
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample genetic
            ancestry group, sex both for adj and raw data.
        - INFO fields for filtering allele frequency (faf) annotations.

    :param subset: Sample subset in dataset. "" is used as a placeholder for the full
        dataset.
    :param description_text: Text describing the sample subset that should be added to
        the INFO description.
    :param data_type: One of "exomes", "genomes", or "joint". Default is "exomes".
    :param pops: Dict of sample global genetic ancestry names for the gnomAD data type.
    :param faf_pops: Dict with gnomAD version (keys) and faf genentic ancestry group
        names (values). Default is FAF_GEN_ANC_GROUPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param label_delimiter: String to use as delimiter when making group label
        combinations. Default is '_'.
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :return: Dictionary containing Subset specific INFO header fields.
    """
    vcf_info_dict = {}
    # Remove unnecessary pop names from FAF_GEN_ANC_GROUPS dict depending on data type
    # and version of FAF_GEN_ANC_GROUPS.
    # Hardcoding faf_pops_version to be 'v1'
    faf_pops_version = 'v1'
    faf_pops_transformed = {pop: GEN_ANC_NAMES[pop] for pop in faf_pops[faf_pops_version]}

    # Add FAF fields to dict.
    faf_label_groups = create_label_groups(
        gen_ancs=list(faf_pops_transformed.keys()),
        sexes=sexes,
        all_groups=['adj'],
    )
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                gen_anc_names=faf_pops_transformed,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                faf=True,
                description_text=description_text,
            ),
        )
    # Add AC, AN, AF, nhomalt fields to dict.
    label_groups = create_label_groups(gen_ancs=list(pops.keys()), sexes=sexes)
    for label_group in label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                gen_anc_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                description_text=description_text,
                callstats=True,
            ),
        )

    # Add grpmax.
    vcf_info_dict.update(
        make_info_dict(
            suffix=subset,
            label_delimiter=label_delimiter,
            gen_anc_names=pops,
            grpmax=True,
            description_text=description_text,
        ),
    )

    # Add fafmax.
    vcf_info_dict.update(
        make_info_dict(
            suffix=subset,
            label_delimiter=label_delimiter,
            gen_anc_names=pops,
            fafmax=True,
            description_text=description_text,
        ),
    )
    if freq_comparison_included:
        ctt_label_groups = create_label_groups(gen_ancs=list(pops.keys()), sexes=sexes)
        for label_group in ctt_label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    freq_ctt=True,
                    label_groups=label_group,
                ),
            )
        vcf_info_dict.update(
            make_info_dict(
                freq_cmh=True,
            ),
        )
        vcf_info_dict.update(
            make_info_dict(
                freq_stat_union=True,
            ),
        )

    return vcf_info_dict


def make_label_combos(
    label_groups: Dict[str, List[str]],
    sort_order: List[str] = SORT_ORDER,
    label_delimiter: str = '_',
) -> List[str]:
    """
    Make combinations of all possible labels for a supplied dictionary of label groups.

    For example, if label_groups is {"sex": ["XY", "XX"], "gen_anc": ["afr", "nfe", "amr"]},
    this function will return ["afr_XY", "afr_XX", "nfe_XY", "nfe_XX", "amr_XY", "amr_XX"]

    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "gen_anc", and value is a list of all possible values for that grouping (e.g. ["XY", "XX"] or ["afr", "nfe", "amr"]).
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :return: list of all possible combinations of values for the supplied label groupings.
    """
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: sort_order.index(x))[0]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x, y in itertools.product(
        anchor_val,
        make_label_combos(copy_label_groups, label_delimiter=label_delimiter),
    ):
        combos.append(f'{x}{label_delimiter}{y}')
    return combos


def make_combo_header_text(
    preposition: str,
    combo_dict: Dict[str, str],
    gen_anc_names: Dict[str, str],
) -> str:
    """
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset.

    For example, if preposition is "for", group_types is ["group", "gen_anc", "sex"], and combo_fields is ["adj", "afr", "XX"],
    this function will return the string " for XX samples in the African-American/African genetic ancestry group".

    :param preposition: Relevant preposition to precede automatically generated text.
    :param combo_dict: Dict with grouping types as keys and values for grouping type as values. This function generates text for these values.
        Possible grouping types are: "group", "gen_anc", "sex", and "subgroup".
        Example input: {"gen_anc": "afr", "sex": "XX"}
    :param gen_anc_names: Dict with global genetic ancestry group names (keys) and genetic ancestry group descriptions (values).
    :return: String with automatically generated description text for a given set of combo fields.
    """
    header_text = ' ' + preposition

    if len(combo_dict) == 1:
        if combo_dict['group'] == 'adj':
            return ''

    if 'sex' in combo_dict:
        header_text = header_text + ' ' + combo_dict['sex']

    header_text = header_text + ' samples'

    if 'subgrp' in combo_dict or 'gen_anc' in combo_dict:
        if 'subgrp' in combo_dict:
            header_text = header_text + f' in the {gen_anc_names[combo_dict["subgrp"]]} genetic ancestry subgroup'

        else:
            header_text = header_text + f' in the {gen_anc_names[combo_dict["gen_anc"]]} genetic ancestry group'

    if 'group' in combo_dict:
        if combo_dict["group"] == 'raw':
            header_text = header_text + ', before removing low-confidence genotypes'

    return header_text


def make_info_dict(
    prefix: str = '',
    suffix: str | None = '',
    prefix_before_metric: bool = True,
    gen_anc_names: Dict[str, str] = GEN_ANC_NAMES,
    label_groups: Dict[str, List[str]] | None = None,
    label_delimiter: str = '_',
    bin_edges: Dict[str, str] | None = None,
    faf: bool = False,
    grpmax: bool = False,
    fafmax: bool = False,
    callstats: bool = False,
    freq_ctt: bool = False,
    freq_cmh: bool = False,
    freq_stat_union: bool = False,
    description_text: str = '',
    age_hist_distribution: str | None = None,
    sort_order: List[str] = SORT_ORDER,
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes of VCF INFO fields.

    Used to populate the INFO fields of the VCF header during export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for grpmax AC, AN, AF, nhomalt, and grpmax genetic ancestry group
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample genetic ancestry group, sex, and subgroup, both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations

    :param prefix: Prefix string for data, e.g. "gnomAD". Default is empty string.
    :param suffix: Suffix string for data, e.g. "gnomAD". Default is empty string.
    :param prefix_before_metric: Whether prefix should be added before the metric (AC, AN, AF, nhomalt, faf95, faf99) in INFO field. Default is True.
    :param gen_anc_names: Dict with global genetic ancestry group names (keys) and genetic ancestry group descriptions (values). Default is GEN_ANC_NAMES.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "gen_anc", and value is a list of all possible values for that grouping (e.g. ["XY", "XX"] or ["afr", "nfe", "amr"]).
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param grpmax: If True, use alternate logic to auto-populate dictionary values associated with grpmax annotations.
    :param fafmax: If True, use alternate logic to auto-populate dictionary values associated with fafmax annotations.
    :param callstats: If True, use alternate logic to auto-populate dictionary values associated with callstats annotations.
    :param freq_ctt: If True, use alternate logic to auto-populate dictionary values associated with frequency contingency table test (CTT) annotations.
    :param freq_cmh: If True, use alternate logic to auto-populate dictionary values associated with frequency Cochran-Mantel-Haenszel (CMH) annotations.
    :param freq_stat_union: If True, use alternate logic to auto-populate dictionary values associated with the union of the contingency table and Cochran-Mantel-Haenszel tests.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :param str age_hist_distribution: Pipe-delimited string of overall age distribution.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    """
    if prefix != '':
        prefix = f'{prefix}{label_delimiter}'
    if suffix != '':
        suffix = f'{label_delimiter}{suffix}'

    info_dict = dict()

    if age_hist_distribution and bin_edges is not None:  # Add bin_edges check:
        age_hist_dict = {
            f'{prefix}age_hist_het_bin_freq{suffix}': {
                'Number': 'A',
                'Description': (
                    f'Histogram of ages of heterozygous individuals{description_text};'
                    f' bin edges are: {bin_edges["het"]}; total number of individuals'
                    f' of any genotype bin: {age_hist_distribution}'
                ),
            },
            f'{prefix}age_hist_het_n_smaller{suffix}': {
                'Number': 'A',
                'Description': (
                    'Count of age values falling below lowest histogram bin edge for'
                    f' heterozygous individuals{description_text}'
                ),
            },
            f'{prefix}age_hist_het_n_larger{suffix}': {
                'Number': 'A',
                'Description': (
                    'Count of age values falling above highest histogram bin edge for'
                    f' heterozygous individuals{description_text}'
                ),
            },
            f'{prefix}age_hist_hom_bin_freq{suffix}': {
                'Number': 'A',
                'Description': (
                    'Histogram of ages of homozygous alternate'
                    f' individuals{description_text}; bin edges are:'
                    f' {bin_edges["hom"]}; total number of individuals of any genotype'
                    f' bin: {age_hist_distribution}'
                ),
            },
            f'{prefix}age_hist_hom_n_smaller{suffix}': {
                'Number': 'A',
                'Description': (
                    'Count of age values falling below lowest histogram bin edge for'
                    f' homozygous alternate individuals{description_text}'
                ),
            },
            f'{prefix}age_hist_hom_n_larger{suffix}': {
                'Number': 'A',
                'Description': (
                    'Count of age values falling above highest histogram bin edge for'
                    f' homozygous alternate individuals{description_text}'
                ),
            },
        }
        info_dict.update(age_hist_dict)

    if grpmax:
        grpmax_dict = {
            f'{prefix}grpmax{suffix}': {
                'Number': 'A',
                'Description': ('Genetic ancestry group with the maximum allele' f' frequency{description_text}'),
            },
            f'{prefix}AC{label_delimiter}grpmax{suffix}': {
                'Number': 'A',
                'Description': (
                    'Allele count in the genetic ancestry group with the maximum allele' f' frequency{description_text}'
                ),
            },
            f'{prefix}AN{label_delimiter}grpmax{suffix}': {
                'Number': 'A',
                'Description': (
                    'Total number of alleles in the genetic ancestry group with the'
                    f' maximum allele frequency{description_text}'
                ),
            },
            f'{prefix}AF{label_delimiter}grpmax{suffix}': {
                'Number': 'A',
                'Description': ('Maximum allele frequency across genetic ancestry' f' groups{description_text}'),
            },
            f'{prefix}nhomalt{label_delimiter}grpmax{suffix}': {
                'Number': 'A',
                'Description': (
                    'Count of homozygous individuals in the genetic ancestry group'
                    f' with the maximum allele frequency{description_text}'
                ),
            },
        }
        info_dict.update(grpmax_dict)

    if fafmax:
        fafmax_dict = {
            f'{prefix}fafmax{label_delimiter}faf95{label_delimiter}max{suffix}': {
                'Number': 'A',
                'Description': (
                    'Maximum filtering allele frequency (using Poisson 95% CI)'
                    f' across genetic ancestry groups{description_text}'
                ),
            },
            f'{prefix}fafmax{label_delimiter}faf95{label_delimiter}max{label_delimiter}gen{label_delimiter}anc{suffix}': {
                'Number': 'A',
                'Description': (
                    'Genetic ancestry group with maximum filtering allele'
                    f' frequency (using Poisson 95% CI){description_text}'
                ),
            },
            f'{prefix}fafmax{label_delimiter}faf99{label_delimiter}max{suffix}': {
                'Number': 'A',
                'Description': (
                    'Maximum filtering allele frequency (using Poisson 99% CI)'
                    f' across genetic ancestry groups{description_text}'
                ),
            },
            f'{prefix}fafmax{label_delimiter}faf99{label_delimiter}max{label_delimiter}gen{label_delimiter}anc{suffix}': {
                'Number': 'A',
                'Description': (
                    'Genetic ancestry group with maximum filtering allele'
                    f' frequency (using Poisson 99% CI){description_text}'
                ),
            },
        }

        info_dict.update(fafmax_dict)

    if callstats or faf or freq_ctt:
        assert label_groups is not None, 'label_groups must be provided if callstats, faf, or freq_ctt is True'
        group_types = sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
        combos = make_label_combos(label_groups, label_delimiter=label_delimiter)

        for combo in combos:
            combo_fields = combo.split(label_delimiter)
            group_dict = dict(zip(group_types, combo_fields))

            for_combo = make_combo_header_text('for', group_dict, gen_anc_names)
            in_combo = make_combo_header_text('in', group_dict, gen_anc_names)

            metrics = ['AC', 'AN', 'AF', 'nhomalt', 'faf95', 'faf99']
            if freq_ctt:
                metrics += ['CTT_odds_ratio', 'CTT_p_value']
            if prefix_before_metric:
                metric_label_dict = {metric: f'{prefix}{metric}{label_delimiter}{combo}{suffix}' for metric in metrics}
            else:
                metric_label_dict = {metric: f'{metric}{label_delimiter}{prefix}{combo}{suffix}' for metric in metrics}

            if callstats:
                combo_dict = {
                    metric_label_dict['AC']: {
                        'Number': 'A',
                        'Description': (f'Alternate allele count{for_combo}{description_text}'),
                    },
                    metric_label_dict['AN']: {
                        'Number': '1',
                        'Description': (f'Total number of alleles{in_combo}{description_text}'),
                    },
                    metric_label_dict['AF']: {
                        'Number': 'A',
                        'Description': (f'Alternate allele frequency{in_combo}{description_text}'),
                    },
                    metric_label_dict['nhomalt']: {
                        'Number': 'A',
                        'Description': ('Count of homozygous' f' individuals{in_combo}{description_text}'),
                    },
                }
            elif faf:
                if ('XX' in combo_fields) | ('XY' in combo_fields):
                    faf_description_text = description_text + ' in non-PAR regions of sex chromosomes only'
                else:
                    faf_description_text = description_text
                combo_dict = {
                    metric_label_dict['faf95']: {
                        'Number': 'A',
                        'Description': (
                            'Filtering allele frequency (using Poisson 95%' f' CI){for_combo}{faf_description_text}'
                        ),
                    },
                    metric_label_dict['faf99']: {
                        'Number': 'A',
                        'Description': (
                            'Filtering allele frequency (using Poisson 99%' f' CI){for_combo}{faf_description_text}'
                        ),
                    },
                }
            else:
                combo_dict = {
                    metric_label_dict['CTT_odds_ratio']: {
                        'Number': 'A',
                        'Description': (
                            "Odds ratio from from Hail's contingency_table_test with"
                            ' `min_cell_count=100` comparing allele frequencies'
                            f' between exomes and genomes{for_combo}{description_text}'
                        ),
                    },
                    metric_label_dict['CTT_p_value']: {
                        'Number': 'A',
                        'Description': (
                            "P-value from Hail's contingency_table_test with"
                            ' `min_cell_count=100` comparing allele frequencies'
                            f' between exomes and genomes{for_combo}{description_text}'
                        ),
                    },
                }
            info_dict.update(combo_dict)
    if freq_cmh:
        cmh_dict = {
            f'{prefix}CMH_chisq{suffix}': {
                'Number': 'A',
                'Description': (
                    'Chi-squared test statistic from the Cochran-Mantel-Haenszel test'
                    ' comparing allele frequencies between exomes and genomes'
                    f' stratified by genetic ancestry group{description_text}'
                ),
            },
            f'{prefix}CMH_p_value{suffix}': {
                'Number': 'A',
                'Description': (
                    'Odds ratio from Cochran-Mantel-Haenszel test comparing allele'
                    ' frequencies between exomes and genomes stratified by genetic'
                    f' ancestry group{description_text}'
                ),
            },
        }
        info_dict.update(cmh_dict)
    if freq_stat_union:
        freq_stat_union_dict = {
            f'{prefix}stat_union_p_value{suffix}': {
                'Number': 'A',
                'Description': (
                    f'p-value from the contingency table or Cochran-Mantel-Haenszel tests{description_text}'
                ),
            },
            f'{prefix}stat_union_test_name{suffix}': {
                'Number': 'A',
                'Description': (
                    f'Name of the test, either contingency_table_test or cochran_mantel_haenszel_test, used to compare allele frequencies between exomes and genomes{description_text}'
                ),
            },
            f'{prefix}stat_union_gen_ancs{suffix}': {
                'Number': '.',
                'Description': (
                    f'List of genetic ancestry groups included in the test. If stat_union_test_name is contingency_table_test, the length of gen_ancs is one and if stat_union_test_name is "cochran_mantel_haenszel_test", the length of "gen_ancs" is greater than one{description_text}'
                ),
            },
        }
        info_dict.update(freq_stat_union_dict)

    return info_dict


def make_hist_bin_edges_expr(
    ht: hl.Table,
    hists: List[str] = HISTS,
    ann_with_hists: Optional[str] = None,
    prefix: str = '',
    label_delimiter: str = '_',
    data_type: str = 'exomes',
    include_age_hists: bool = True,
    for_joint: bool = False,
) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :param ann_with_hists: Name of row annotation containing histogram data. In exomes or
        genomes release HT, `histograms` is a row, but in the joint release HT, it's
        under the row of `exomes`, `genomes`, or `joint`.
    :param prefix: Prefix text for age histogram bin edges.  Default is empty string.
    :param label_delimiter: String used as delimiter between prefix and histogram annotation.
    :param include_age_hists: Include age histogram annotations.
    :return: Dictionary keyed by histogram annotation name, with corresponding
        reformatted bin edges for values.
    """
    # Add underscore to prefix if it isn't empty
    if prefix:
        prefix += label_delimiter

    edges_dict = {}

    # NOTE: `age_hists` is not a struct in our exomes and genomes tables. `age_hist_het` and `age_hist_hom` are
    # separate top-level structs not under `histograms`
    if include_age_hists:
        for call_type in ['het', 'hom']:
            if ann_with_hists:
                bin_edges = (
                    ht.filter(
                        hl.is_defined(
                            ht[ann_with_hists].histograms.age_hists[f'age_hist_{call_type}'].bin_edges,
                        ),
                    )[ann_with_hists]
                    .histograms.age_hists[f'age_hist_{call_type}']
                    .bin_edges.take(1)[0]
                )
            else:
                if for_joint:
                    bin_edges = (
                        ht.filter(
                            hl.is_defined(
                                ht.histograms.age_hists[f'age_hist_{call_type}'].bin_edges,
                            ),
                        )
                        .histograms.age_hists[f'age_hist_{call_type}']
                        .bin_edges.take(1)[0]
                    )
                else:
                    bin_edges = ht.filter(
                        hl.is_defined(
                            ht[f'age_hist_{call_type}'].bin_edges,
                        ),
                    )[
                        f'age_hist_{call_type}'
                    ].bin_edges.take(1)[0]

            if bin_edges:
                edges_dict[f'{prefix}{call_type}'] = '|'.join(
                    map(lambda x: f'{x:.1f}', bin_edges),
                )

    for hist in hists:
        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in [f'{prefix}raw_qual_hists', f'{prefix}qual_hists']:
            hist_name = hist if 'raw' not in hist_type else f'{prefix}{hist}_raw'

            # This if-statement isn't going to matter because neither of our exomes
            # or genomes have histograms stored under ht.exomes/genomes/joint
            # But leaving it in for future-proofing.
            if ann_with_hists:
                bin_edges = (
                    ht.filter(
                        hl.is_defined(
                            ht[ann_with_hists].histograms[hist_type][hist].bin_edges,
                        ),
                    )[ann_with_hists]
                    .histograms[hist_type][hist]
                    .bin_edges.take(1)[0]
                )
            else:
                bin_edges = (
                    ht.filter(hl.is_defined(ht.histograms[hist_type][hist].bin_edges))
                    .histograms[hist_type][hist]
                    .bin_edges.take(1)[0]
                )
            if bin_edges:
                edges_dict[hist_name] = '|'.join(
                    map(
                        lambda x: f'{x:.2f}' if 'ab' in hist else str(int(x)),
                        bin_edges,
                    ),
                )

    return edges_dict


def make_hist_dict(
    bin_edges: Dict[str, str],
    adj: bool,
    hist_metric_list: List[str] = HISTS,
    label_delimiter: str = '_',
    drop_n_smaller_larger: bool = False,
    prefix: str = '',
    suffix: str = '',
    description_text: str = '',
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations.

    :param bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values.
    :param adj: Whether to create a header dict for raw or adj quality histograms.
    :param hist_metric_list: List of hists for which to build hist info dict
    :param label_delimiter: String used as delimiter in values stored in hist_metric_list.
    :param drop_n_smaller_larger: Whether to drop n_smaller and n_larger annotations from header dict. Default is False.
    :param prefix: Prefix text for histogram annotations. Default is empty string.
    :param suffix: Suffix text for histogram annotations. Default is empty string.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes.
    """
    if prefix != '':
        prefix = f'{prefix}{label_delimiter}'
    if suffix != '':
        suffix = f'{label_delimiter}{suffix}'

    header_hist_dict = {}
    for hist in hist_metric_list:
        # Get hists for both raw and adj data
        # Add '_raw' to quality histograms calculated on raw data
        if not adj:
            hist = f'{hist}_raw'

        edges = bin_edges[hist]
        hist_fields = hist.split(label_delimiter)
        hist_text = hist_fields[0].upper()

        if hist_fields[2] == 'alt':
            hist_text = hist_text + ' in heterozygous individuals'
        if adj:
            hist_text = hist_text + ' calculated on high quality genotypes'

        hist_dict = {
            f'{prefix}{hist}_bin_freq{suffix}': {
                'Number': 'A',
                'Description': (f'Histogram for {hist_text}{description_text}; bin edges are:' f' {edges}'),
            },
        }
        # These annotations are frequently zero and are dropped from gnomad
        # releases for most histograms.
        if not drop_n_smaller_larger:
            hist_dict.update(
                {
                    f'{prefix}{hist}_n_smaller{suffix}': {
                        'Number': 'A',
                        'Description': (
                            f'Count of {hist_fields[0].upper()} values falling below'
                            f' lowest histogram bin edge {hist_text}{description_text}'
                        ),
                    },
                    f'{prefix}{hist}_n_larger{suffix}': {
                        'Number': 'A',
                        'Description': (
                            f'Count of {hist_fields[0].upper()} values falling above'
                            f' highest histogram bin edge {hist_text}{description_text}'
                        ),
                    },
                },
            )
        # Only add n_larger for dp qual histograms.
        if 'dp' in hist:
            hist_dict.update(
                {
                    f'{prefix}{hist}_n_larger{suffix}': {
                        'Number': 'A',
                        'Description': (
                            f'Count of {hist_fields[0].upper()} values falling above'
                            f' highest histogram bin edge {hist_text}{description_text}'
                        ),
                    },
                },
            )

        header_hist_dict.update(hist_dict)

    return header_hist_dict


def add_as_info_dict(
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT,
    as_fields: List[str] = AS_FIELDS,
) -> Dict[str, Dict[str, str]]:
    """
    Update info dictionary with allele-specific terms and their descriptions.

    Used in VCF export.

    :param info_dict: Dictionary containing site-level annotations and their descriptions. Default is INFO_DICT.
    :param as_fields: List containing allele-specific fields to be added to info_dict. Default is AS_FIELDS.
    :return: Dictionary with allele specific annotations, their descriptions, and their VCF number field.
    """
    as_dict: Dict[str, Dict[str, str]] = {}
    for field in as_fields:
        try:
            # Strip AS_ from field name
            site_field = field[3:]

            # Get site description from info dictionary and make first letter lower case
            first_letter = info_dict[site_field]['Description'][0].lower()
            rest_of_description = info_dict[site_field]['Description'][1:]

            as_dict[field] = {}
            as_dict[field]['Number'] = 'A'
            as_dict[field]['Description'] = f'Allele-specific {first_letter}{rest_of_description}'

        except KeyError:
            logger.warning(f'{field} is not present in input info dictionary!')

    return as_dict


def populate_info_dict(
    info_fields: List[str],
    bin_edges: Dict[str, str],
    age_hist_distribution: str | None = None,
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT,
    subset_list: List[str] = SUBSETS['exomes'],
    pops: Dict[str, str] = GEN_ANC_GROUPS['exomes'],  # type: ignore[assignment]
    faf_pops: Dict[str, List[str]] = FAF_GEN_ANC_GROUPS,
    sexes: List[str] = SEXES,
    in_silico_dict: Dict[str, Dict[str, str]] | None = IN_SILICO_ANNOTATIONS_INFO_DICT,
    vrs_fields_dict: Dict[str, Dict[str, str]] | None = VRS_FIELDS_DICT,
    label_delimiter: str = '_',
    data_type: str = 'exomes',
    freq_comparison_included: bool = False,
    extra_suffix: str | None = None,
    extra_description_text: str | None = None,
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` and `make_hist_dict` to populate INFO dictionary.

    Used during VCF export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for
          heterozygous and homozygous variant carriers).
        - INFO fields for grpmax AC, AN, AF, nhomalt, and grpmax genetic ancestry group.
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample genetic
          ancestry group, sex both for adj and raw data.
        - INFO fields for filtering allele frequency (faf) annotations.
        - INFO fields for variant histograms (hist_bin_freq for each histogram and
          hist_n_larger for DP histograms).

    :param info_fields: List of info fields to add to the info dict. Default is None.
    :param bin_edges: Dictionary of variant annotation histograms and their associated
        bin edges.
    :param age_hist_distribution: Pipe-delimited string of overall age histogram bin
        frequency.
    :param info_dict: INFO dict to be populated.
    :param subset_list: List of sample subsets in dataset. Default is SUBSETS["exomes"].
    :param pops: Dict of sample global genetic ancestry names for the gnomAD data type.
    :param faf_pops: Dict with gnomAD version (keys) and faf genentic ancestry group
        names (values). Default is FAF_GEN_ANC_GROUPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param in_silico_dict: Dictionary of in silico predictor score descriptions.
    :param vrs_fields_dict: Dictionary with VRS annotations.
    :param label_delimiter: String to use as delimiter when making group label
        combinations.
    :param data_type: Data type to populate info dict for. One of "exomes" or
        "genomes". Default is "exomes".
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :param extra_suffix: Suffix to add to INFO field.
    :param extra_description_text: Extra description text to add to INFO field.
    :return: Updated INFO dictionary for VCF export.
    """
    vcf_info_dict = {}
    if data_type == 'joint':
        # vcf_info_dict stays empty if data_type is 'joint' and subset is not 'joint'
        if 'joint' in subset_list:
            vcf_info_dict.update(JOINT_REGION_FLAGS_INFO_DICT)
    else:
        # Get existing info fields from predefined info_dict, e.g. `FS`,
        # `non_par`, `negative_train_site`...
        vcf_info_dict.update(info_dict)
        vcf_info_dict = {f: vcf_info_dict[f] for f in info_fields if f in vcf_info_dict}
        # Add allele-specific fields to info dict, including AS_VQSR_FIELDS
        vcf_info_dict.update(
            add_as_info_dict(info_dict=info_dict, as_fields=AS_FIELDS + AS_VQSR_FIELDS),
        )

    for subset in subset_list:
        subset_pops = deepcopy(pops)
        if data_type == 'joint':
            description_text = f' in {subset} dataset' if subset != '' else ''
        else:
            description_text = '' if subset == '' else f' in {subset} subset'

        vcf_info_dict.update(
            populate_subset_info_dict(
                subset=subset,
                description_text=description_text,
                data_type=data_type,
                pops=subset_pops,
                faf_pops=faf_pops,
                sexes=sexes,
                label_delimiter=label_delimiter,
                freq_comparison_included=freq_comparison_included,
            ),
        )

    if age_hist_distribution:
        age_hist_distribution = '|'.join(str(x) for x in age_hist_distribution)

    # Add age histogram data to info dict.
    vcf_info_dict.update(
        make_info_dict(
            suffix=extra_suffix if extra_suffix else '',
            label_delimiter=label_delimiter,
            bin_edges=bin_edges,
            age_hist_distribution=age_hist_distribution,
            description_text=extra_description_text if extra_description_text else '',
        ),
    )

    # Add variant quality histograms to info dict.
    vcf_info_dict.update(
        make_hist_dict(
            bin_edges,
            adj=True,
            drop_n_smaller_larger=True,
            suffix=extra_suffix if extra_suffix else '',
            description_text=extra_description_text if extra_description_text else '',
        ),
    )
    if data_type == 'joint':
        vcf_info_dict.update(JOINT_FILTERS_INFO_DICT)
        return vcf_info_dict

    return vcf_info_dict


def adjust_interval_padding(ht: hl.Table, padding: int) -> hl.Table:
    """
    Adjust interval padding in HT.

    .. warning::

        This function can lead to overlapping intervals, so it is not recommended for
        most applications. For example, it can be used to filter a variant list to all
        variants within the returned interval list, but would not work for getting an
        aggregate statistic for each interval if the desired output is independent
        statistics.

    :param ht: HT to adjust.
    :param padding: Padding to use.
    :return: HT with adjusted interval padding.
    """
    return ht.key_by(
        interval=hl.locus_interval(
            ht.interval.start.contig,
            ht.interval.start.position - padding,
            ht.interval.end.position + padding,
            # Include the end of the intervals to capture all variants.
            reference_genome=ht.interval.start.dtype.reference_genome,
            includes_end=True,
        ),
    )


def add_capture_region_flags(ht: hl.Table) -> hl.Table:
    capture_intervals_paths: list[str] = config_retrieve(['large_cohort', 'browser', 'capture_intervals'])
    capture_tables = [
        hl.import_locus_intervals(str(path), reference_genome=genome_build()) for path in capture_intervals_paths
    ]
    if len(capture_intervals_paths) > 1:
        # Deduplicate keys, keeping exactly one row for each unique key (key='interval').
        capture_intervals_ht = hl.Table.union(*capture_tables)
    else:
        capture_intervals_ht = capture_tables[0]

    calling_intervals_path = config_retrieve(['large_cohort', 'browser', 'calling_intervals'])
    calling_interval_padding = config_retrieve(['large_cohort', 'browser', 'calling_interval_padding'], default=150)
    calling_intervals_ht = hl.import_locus_intervals(str(calling_intervals_path), reference_genome=genome_build())

    intervals: list[tuple[str, hl.Table]] = [
        ('calling', adjust_interval_padding(calling_intervals_ht, calling_interval_padding)),
    ] + [
        ('capture', capture_intervals_ht),
    ]

    new_flags = {
        **{f"outside_{c}_region": hl.is_missing(i[ht.locus]) for c, i in intervals},
    }

    ht = ht.annotate(region_flags=ht.region_flags.annotate(**new_flags))
    return ht


def repartition_frequencies_table(
    ht_path: str,
    data_type: str,
):
    """
    Repartition the genomes frequency table.
    # TODO fix this readme
    """
    if data_type == "exome":
        return(ht_path)
    else:
        repartitioned_path = output_path(f"{data_type}_repartitioned.ht", category="tmp")
        freq_table = hl.read_table(ht_path).repartition(n = 10000).checkpoint(repartitioned_path, overwrite=True)
        return(repartitioned_path)


def run_browser_vcf_data_download(
    ht_path: hl.Table,
    data_type: str,
    contig: str,
    vcf_outpath: str,
    joint_included: bool = False,
    vqsr_ht_path: str | None = None,
    exome_freq_ht_path: str | None = None,
    genome_freq_ht_path: str | None = None,
):
    """
    Export gnomAD frequency data to VCF format with appropriate header information.

    NOTE: For joint data export, uses the processed joint table from browser prepare
    rather than the raw joint frequencies table. For exomes/genomes, uses the
    original frequencies stage output.
    """
    ht = hl.read_table(ht_path)

    if vqsr_ht_path and data_type != 'joint':
        # genomes and exomes frequency tables have allele_info dropped
        # Annotating back allele_info to have allele_info fields in info dict
        vqsr_ht = hl.read_table(vqsr_ht_path)
        ht = ht.annotate(allele_info=vqsr_ht[ht.key].allele_info)

    if data_type == 'exomes':
        # annotate capture and calling regions to ht.region_flags
        ht = add_capture_region_flags(ht)

    score_cutoffs = {
        seq_type: process_score_cutoffs(
            snv_bin_cutoff=config_retrieve(['large_cohort', 'browser', f'snp_bin_threshold_{seq_type}']),
            indel_bin_cutoff=config_retrieve(['large_cohort', 'browser', f'indel_bin_threshold_{seq_type}']),
            aggregated_bin_ht_path=config_retrieve(
                ['large_cohort', 'browser', f'aggregated_bin_ht_path_{seq_type}'],
            ),
            snv_bin_id=config_retrieve(['large_cohort', 'browser', f'snp_bin_id_{seq_type}'], 'bin'),
            indel_bin_id=config_retrieve(['large_cohort', 'browser', f'indel_bin_id_{seq_type}'], 'bin'),
        )
        for seq_type in ['exome', 'genome']
    }

    for_joint = data_type == 'joint'

    validate_hts = {}

    if data_type == 'joint':
        iter_data_types = ['exomes', 'genomes', 'joint']
    else:
        iter_data_types = [data_type]

    # Pre-add filter fields for joint.exomes.filters and joint.genomes.filters from
    # exomes_freq_ht and genomes_freq_ht respectively. Except for AC0
    if for_joint and (exome_freq_ht_path is not None) and (genome_freq_ht_path is not None):
        exomes_freq_ht = hl.read_table(exome_freq_ht_path)
        genomes_freq_ht = hl.read_table(genome_freq_ht_path)
        exomes_freq_ht = get_filters_expr(ht=exomes_freq_ht, score_cutoffs=score_cutoffs['exome'])
        genomes_freq_ht = get_filters_expr(ht=genomes_freq_ht, score_cutoffs=score_cutoffs['genome'])

        # Annotate back the filters to the joint.exomes.filters and joint.genomes.filters
        ht = ht.annotate(
            exomes=ht.exomes.annotate(filters=exomes_freq_ht[ht.key].filters),
            genomes=ht.genomes.annotate(filters=genomes_freq_ht[ht.key].filters),
        )

    for dt in iter_data_types:
        if for_joint:
            dt_ht = select_type_from_joint_ht(ht, dt)
        else:
            dt_ht = ht

        logger.info(f'Preparing {dt} HT for validity checks and export...')
        # For joint validation: annotate filters from original exomes/genomes tables since
        # joint HT lacks filters after subsetting by select_type_from_joint_ht()
        dt_ht, rename_dict = prepare_ht_for_validation(
            dt_ht,
            data_type=dt,
            joint_included=joint_included,
            freq_comparison_included=(dt == 'joint'),
            for_joint_validation=for_joint,
            score_cutoffs=score_cutoffs[data_type[:-1]] if not for_joint else None,  # de-pluralise
        )

        if for_joint:
            ordered_rename_dict = {key: rename_dict.get(key, key) for key in dt_ht.info.keys()}
            dt_ht = dt_ht.annotate(info=dt_ht.info.rename(ordered_rename_dict))
            if dt != 'joint':
                dt_ht = dt_ht.annotate(
                    info=dt_ht.info.annotate(**{f'{dt}_filters': dt_ht.filters}),
                )
                dt_ht = dt_ht.select('info')

            dt_ht = dt_ht.select_globals(
                **{f'{dt}_{f}': dt_ht[f] for f in dt_ht.globals},
            )
        validate_hts[dt] = dt_ht

    ht = validate_hts[data_type]
    if for_joint:
        in_joint_ht = set(ht.info.keys())
        for dt in ['exomes', 'genomes']:
            info_expr = validate_hts[dt][ht.key].info
            info_expr = info_expr.select(
                *[f for f in info_expr if f not in in_joint_ht],
            )
            ht = ht.annotate(info=ht.info.annotate(**info_expr))
            ht = ht.annotate_globals(**validate_hts[dt].index_globals())

        ht = get_joint_filters(ht)

    logger.info(f'Checkpointing validated_ht to {output_path(f"{contig}_{data_type}_validated.ht", category="tmp")}...')
    validated_ht = ht.checkpoint(output_path(f"{contig}_{data_type}_validated.ht", category="tmp"), overwrite=True)

    ht = hl.read_table(ht_path)
    # Our VQSR features match gnomAD's filtering_model structure, but cutoff values
    # aren't available in frequency tables - they're added during browser prepare. So adding them here.
    # Only applicable to genomes and exomes, not joint
    if data_type != 'joint':
        ht = ht.annotate_globals(
            filtering_model=hl.struct(
                filter_name='AS_VQSR',
                score_name='AS_VQSLOD',
                snv_cutoff=hl.struct(
                    bin=score_cutoffs[data_type[:-1]]['snv']['min_score'],
                    min_score=score_cutoffs[data_type[:-1]]['snv']['min_score'],
                ),
                indel_cutoff=hl.struct(
                    bin=score_cutoffs[data_type[:-1]]['indel']['min_score'],
                    min_score=score_cutoffs[data_type[:-1]]['indel']['min_score'],
                ),
                snv_training_variables=[
                    'AS_QD',
                    'AS_MQRankSum',
                    'AS_ReadPosRankSum',
                    'AS_FS',
                    'AS_SOR',
                    'AS_MQ',
                ],
                indel_training_variables=[
                    'AS_QD',
                    'AS_MQRankSum',
                    'AS_ReadPosRankSum',
                    'AS_FS',
                    'AS_SOR',
                ],
            ),
        )
        inbreeding_coeff_cutoff = hl.float64(
            config_retrieve(['large_cohort', 'browser', 'inbreeding_coeff_cutoff']),
        )
        ht = ht.annotate_globals(
            inbreeding_coeff_cutoff=inbreeding_coeff_cutoff,
        )

    if not for_joint:
        # v4 Genomes drops subsets from VCF
        subsets = SUBSETS['exomes'] if data_type == 'exomes' else []
        header_dict = prepare_vcf_header_dict(
            ht,
            validated_ht=validated_ht,
            info_fields=list(validated_ht.info),
            bin_edges=make_hist_bin_edges_expr(
                ht,
                data_type=data_type,
                include_age_hists=True,
            ),
            age_hist_distribution=hl.eval(ht.age_distribution.bin_freq),
            subset_list=subsets,
            pops=GEN_ANC_GROUPS[data_type],  # type: ignore[arg-type]
            data_type=data_type,
            joint_included=joint_included,
        )
    else:
        header_dict = {'filter': make_vcf_filter_dict(joint=True), 'info': {}}  # type: ignore[dict-item]
        for dt in ['exomes', 'genomes', 'joint']:
            dt_ht = select_type_from_joint_ht(ht, dt)
            temp_header_dict = prepare_vcf_header_dict(
                dt_ht,
                validated_ht=validated_ht,
                info_fields=[f for f in validated_ht.info.keys() if dt in f],
                bin_edges=make_hist_bin_edges_expr(
                    dt_ht,
                    data_type=data_type,
                    include_age_hists=True,
                    for_joint=for_joint,
                ),
                age_hist_distribution=hl.eval(dt_ht.age_distribution.bin_freq),
                subset_list=[dt],
                pops=GEN_ANC_GROUPS[dt],  # type: ignore[arg-type]
                data_type='joint',
                joint_included=joint_included,
                freq_comparison_included=(dt == 'joint'),
                extra_suffix=dt,
                extra_description_text=f' in {dt} dataset',
            )
            header_dict['info'].update(temp_header_dict)  # type: ignore[arg-type]

    with hl.hadoop_open(output_path(f'{contig}_{data_type}_header_dict.pkl', category='tmp'), 'wb') as p:
        pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

    ht = validated_ht

    if contig:
        ht = hl.filter_intervals(
            ht,
            [hl.parse_locus_interval(contig, reference_genome='GRCh38')],
        )

    ht, new_row_annots = format_validated_ht_for_export(ht, data_type=data_type)

    ordered_vcf_info_dict = {f: header_dict['info'][f] for f in list(ht.info) if f in header_dict['info']}
    header_dict.update({'info': ordered_vcf_info_dict})

    hl.export_vcf(
        ht,
        vcf_outpath,
        metadata=header_dict,
        # append_to_header=append_to_vcf_header_path(data_type=data_type),
        tabix=True,
    )
