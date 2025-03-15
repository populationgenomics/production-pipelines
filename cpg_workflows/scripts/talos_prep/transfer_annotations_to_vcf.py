#!/usr/bin/env python3

"""
Takes the full un-annotated joint-callset VCF,
A HailTable of the formatted annotations,
Integrates the two, writing a MatrixTable representation of the fully annotated VCF
"""

import logging
from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import init_batch


def cli_main():
    """
    take an input VCF and an output MT path
    also supply the alpha_missense table created by parse_amissense_into_ht.py
    """

    parser = ArgumentParser(description='Takes a HT of annotations, and a callset VCF, and combines into a MT')
    parser.add_argument('--input_path', help='Path to the input, MatrixTable or VCF', required=True)
    parser.add_argument('--annotations', help='Path to the annotation HT', required=True)
    parser.add_argument('--output', help='output Table path, must have a ".mt" extension', required=True)
    args = parser.parse_args()

    assert args.output.endswith('.mt'), 'Output path must end in .mt'

    main(
        input_path=args.input_path,
        output_path=args.output,
        annotations=args.annotations,
    )


def main(
    input_path: str,
    output_path: str,
    annotations: str,
):
    """
    Takes a Hail-Table of annotations, a joint-called VCF, reads the VCF as a MatrixTable and hops the annotations over

    Args:
        input_path (str): path to the full callset VCF or MatrixTable
        output_path (str): path to write the resulting MatrixTable to
        annotations (str): path to a Hail Table containing annotations
    """
    init_batch()
    # hl.context.init_spark(master='local[8]', default_reference='GRCh38')

    if input_path.endswith('.mt'):
        # read the MatrixTable directly
        mt = hl.read_matrix_table(input_path)
    elif input_path.endswith(('.vcf.gz', '.vcf.bgz')):
        # read the VCF into a MatrixTable
        mt = hl.import_vcf(
            input_path,
            array_elements_required=False,
            force_bgz=True,
        )
    else:
        raise ValueError('Input path must end in .mt, .vcf.gz, or .vcf.bgz')

    # read the annotations into a Table
    ht = hl.read_table(annotations)

    # syntax sweeter for later on
    matched_annotations = ht[mt.row_key]

    # take the fuller annotations in gnomAD 4.1 joint exomes ~ genomes
    mt = mt.annotate_rows(
        gnomad=hl.struct(
            gnomad_AC=matched_annotations.info.gnomad_AC_joint,
            gnomad_AF=matched_annotations.info.gnomad_AF_joint,
            gnomad_AN=matched_annotations.info.gnomad_AN_joint,
            gnomad_AC_XY=matched_annotations.info.gnomad_AC_joint_XY,
            gnomad_AF_XY=matched_annotations.info.gnomad_AF_joint_XY,
            gnomad_FAF=matched_annotations.info.gnomad_faf_95_joint,
            gnomad_HomAlt=matched_annotations.info.gnomad_HomAlt_joint,
        ),
        transcript_consequences=matched_annotations.transcript_consequences,
        gene_ids=matched_annotations.gene_ids,
    )

    mt.describe()

    mt.write(output_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    cli_main()
