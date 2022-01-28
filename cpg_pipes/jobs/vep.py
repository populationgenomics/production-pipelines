#!/usr/bin/env python3

"""
Batch pipeline to laod data into seqr
"""
from typing import Optional

from hailtop.batch.job import Job

from cpg_pipes import resources, hailbatch


def vep(
    b,
    vcf_path: str,
    out_vcf_path: Optional[str] = None,
) -> Job:
    """
    Runs VEP on provided VCF.
    """
    j = b.new_job('VEP')
    j.image(resources.VEP_IMAGE)
    hailbatch.STANDARD.set_resources(j, storage_gb=25, mem_gb=20)

    cmd = f"""\
    LOFTEE_DIR=/io/batch/loftee
    CACHE_DIR=/io/batch/cache
    mkdir -p $LOFTEE_DIR
    mkdir -p $CACHE_DIR
    gsutil cat {resources.VEP_LOFTEE} | tar -xf - -C $LOFTEE_DIR/
    ls $LOFTEE_DIR
    du -sh $LOFTEE_DIR

    gsutil cat {resources.VEP_CACHE} | tar -xf - -C $CACHE_DIR/
    ls $CACHE_DIR
    du -sh $CACHE_DIR
    FASTA=$CACHE_DIR/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    ls $(dirname $FASTA)

    gsutil cp {vcf_path} input.vcf.gz
    
    GEPR_PATH=$LOFTEE_DIR/gerp_conservation_scores.homo_sapiens.GRCh38.bw
    GEPR_PATH=

    vep \\
    --vcf \\
    --format vcf \\
    -i input.vcf.gz \\
    --everything \\
    --allele_number \\
    --no_stats \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache $CACHE_DIR/vep/ \\
    --dir_plugins /root/micromamba/share/ensembl-vep-105.0-0/ \\
    --fasta $FASTA \\
    --plugin LoF,loftee_path:/root/micromamba/share/ensembl-vep-105.0-0/,gerp_bigwig:$LOFTEE_DIR/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:$LOFTEE_DIR/human_ancestor.fa.gz,conservation_file:$LOFTEE_DIR/loftee.sql \\
    -o {j.out_vcf} \\
    --compress_output bgzip
    """
    j.command(hailbatch.wrap_command(
        cmd, 
        setup_gcp=True, 
        monitor_space=True, 
        define_retry_function=True
    ))
    if out_vcf_path:
        b.write_output(j.out_vcf, out_vcf_path)
    return j
