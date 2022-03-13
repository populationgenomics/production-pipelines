"""
Hail Batch Germline CNV caller based on:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531152
"""

from typing import List

import hailtop.batch as hb


def main(reference: str, crams: List[str]):

    batch = hb.Batch("gvnc")

    _reference = batch.read_input_group(
        base=reference,
        fai=reference + ".fai",
        dict=reference.replace(".fasta", "").replace(".fna", "").replace(".fa", "")
        + ".dict",
    )
    _crams = [batch.read_input_group(base=cram, crai=cram + ".crai") for cram in crams]

    interval_list = preprocess_intervals(batch, _reference)

    cram_counts = []
    for cram in _crams:
        cram_counts.append(collect_read_counts(batch, cram))

    ploidy_priors = determine_germline_contig_ploidy(batch, interval_list, cram_counts)
    output_vcf = run_germline_cnv_caller(
        batch, interval_list, cram_counts, ploidy_priors
    )


def preprocess_intervals(batch, reference):
    j = batch.new_job("preprocess-intervals")
    j.command(
        f"""\
gatk PreprocessIntervals \
        -R {reference.fasta} \
        --padding 0 \
        -imr OVERLAPPING_ONLY \
        -O {j.out_interval_list}
"""
    )

    return j.out_interval_list


def collect_read_counts(batch, cram):

    j = batch.new_job("collect_read_count_cram")

    out_counts_tsv = j.out_counts_tsv

    j.command(
        f"""
gatk CollectReadCounts \
    -L chr20sub.interval_list \
    -R ref/Homo_sapiens_assembly38.fasta \
    -imr OVERLAPPING_ONLY \
    -I {cram.cram} \
    --format TSV \
    -O {out_counts_tsv}
    """
    )

    return out_counts_tsv


def determine_germline_contig_ploidy(batch, interval_list, counts: List):
    j = batch.new_job("determine-germline-contig-ploidy")

    counts_str = "\\\n\t".join(f"-I {count}" for count in counts)

    # declare outputs
    out_ploidy_priors = j.out_ploidy_priors

    # command
    j.command(
        f"""
gatk DetermineGermlineContigPloidy \
    -L {interval_list} \
    --interval-merging-rule OVERLAPPING_ONLY \
    {counts_str}
    --contig-ploidy-priors {out_ploidy_priors} \
    --output . \
    --output-prefix ploidy \
    --verbosity DEBUG
"""
    )

    # this produces a LOT more files than originally expected

    return out_ploidy_priors


def run_germline_cnv_caller(
    batch, interval_list, counts, ploidy_priors, annotated_intervals
):

    j = batch.new_job("germline-cnv-caller")

    counts_str = "\\\n\t".join(f"-I {count}" for count in counts)

    j.command(
        f"""
gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L {interval_list} \
    {counts_str} \
    --contig-ploidy-calls ploidy-calls \
    --annotated-intervals twelveregions.annotated.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output out \
    --output-prefix na- \
    --verbosity DEBUG"""
    )

    return
