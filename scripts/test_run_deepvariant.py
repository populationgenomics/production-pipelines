#!/usr/bin/env python3

from typing import TYPE_CHECKING

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def run(output_path: str) -> 'Job':
    """
    This is a simple example of a job that writes a statement to a file.

    Args:
        statement (str): the intended file contents
        output_file (str): the path to write the file to

    Returns:
        the resulting job
    """

    # create a job
    j = get_batch().new_job('DeepVariant')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(image_path('deepvariant'))

    # copy test data
    j.command(
        f"""
        INPUT_DIR="$BATCH_TMPDIR/quickstart-testdata"
        CHECKPOINT_DIR="$BATCH_TMPDIR/checkpoint"
        DATA_HTTP_DIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"

        mkdir -p ${{INPUT_DIR}}
        mkdir -p ${{CHECKPOINT_DIR}}
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/NA12878_S1.chr20.10_10p1mb.bam"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/NA12878_S1.chr20.10_10p1mb.bam.bai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.bed"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.fai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz.fai"
        wget -P ${{INPUT_DIR}} "${{DATA_HTTP_DIR}}/ucsc.hg19.chr20.unittest.fasta.gz.gzi"

        # Run make_examples
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --vcf_stats_report=true \
        --ref=${{INPUT_DIR}}/ucsc.hg19.chr20.unittest.fasta \
        --reads=${{INPUT_DIR}}/NA12878_S1.chr20.10_10p1mb.bam \
        --regions "chr20:10,000,000-10,010,000" \
        --output_vcf=/output/output.vcf.gz \
        --output_gvcf={j.outfile} \
        --intermediate_results_dir /output/intermediate_results_dir \
        --num_shards=1
        """,
    )

    # write the output to the expected location
    get_batch().write_output(j.outfile, output_path)

    # return the job
    return j


def run_pangenome_aware(output_vcf: str, output_gvcf: str) -> 'Job':
    # create a job
    j = get_batch().new_job('DeepVariant PangenomeAware')

    # choose an image to run this job in (default is bare ubuntu)
    j.image(image_path('deepvariant_pangenome_aware'))
    j.memory('8Gi')

    # copy test data
    j.command(
        f"""
        set -ex
        mkdir -p reference
        FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
        curl ${{FTPDIR}}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gunzip > reference/GRCh38_no_alt_analysis_set.fasta
        curl ${{FTPDIR}}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > reference/GRCh38_no_alt_analysis_set.fasta.fai
        ls -l reference

        mkdir -p benchmark
        FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38
        curl ${{FTPDIR}}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
        curl ${{FTPDIR}}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
        curl ${{FTPDIR}}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi > benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
        ls -l benchmark

        mkdir -p input
        HTTPDIR=https://storage.googleapis.com/deepvariant/case-study-testdata
        curl ${{HTTPDIR}}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam
        curl ${{HTTPDIR}}/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai > input/HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam.bai
        ls -l input

        mkdir -p input
        ls -l input
        HTTPDIR=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38
        curl -L "${{HTTPDIR}}/hprc-v1.1-mc-grch38.gbz" -o input/hprc-v1.1-mc-grch38.gbz

        # Run pangenome-aware DeepVariant
        mkdir -p output
        mkdir -p output/intermediate_results_dir
        /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
        --model_type WGS \
        --ref /reference/GRCh38_no_alt_analysis_set.fasta \
        --reads /input//HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam  \
        --pangenome /input/hprc-v1.1-mc-grch38.gbz \
        --output_vcf {j.out_vcf} \
        --output_gvcf {j.out_gvcf} \
        --num_shards $(nproc) \
        --regions chr20 \
        --intermediate_results_dir /output/intermediate_results_dir
        """,
    )

    # write the output to the expected location
    get_batch().write_output(j.out_vcf, output_vcf)
    get_batch().write_output(j.out_gvcf, output_gvcf)

    # return the job
    return j


def main():
    """
    Main entry point for testing the run function.
    """
    init_batch()
    output_vcf = 'gs://cpg-bioheart-test/deepvariant_test/test_out.vcf.gz'
    output_gvcf = 'gs://cpg-bioheart-test/deepvariant_test/test_out.g.vcf.gz'
    job = run_pangenome_aware(output_vcf, output_gvcf)
    print(f"Job {job.name} created with output paths: {output_vcf}, {output_gvcf}")

    get_batch().run()


if __name__ == '__main__':
    main()
