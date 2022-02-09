#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""

from typing import Optional

from hailtop.batch.job import Job

from cpg_pipes import images, ref_data, files
from cpg_pipes.hb.command import wrap_command


def vep(
    b,
    vcf_path: str,
    out_vcf_path: Optional[str] = None,
    overwrite: bool = True
) -> Job:
    """
    Runs VEP on provided VCF.
    """
    j = b.new_job('VEP')
    if out_vcf_path and files.can_reuse(out_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.VEP_IMAGE)
    j.storage('50G')
    j.memory('50G')

    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': '$LOFTEE_DIR/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': '$LOFTEE_DIR/human_ancestor.fa.gz',
        'conservation_file': '$LOFTEE_DIR/loftee.sql',
    }
    
    cmd = f"""\
    CACHE_DIR=/io/batch/cache
    LOFTEE_DIR=/io/batch/loftee
    LOFTEE_PLUGIN_PATH=/root/micromamba/share/ensembl-vep-105.0-0
    FASTA=$CACHE_DIR/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    mkdir -p $CACHE_DIR
    mkdir -p $LOFTEE_DIR
    gsutil cat {ref_data.VEP_CACHE} | tar -xf - -C $CACHE_DIR/
    gsutil cat {ref_data.VEP_LOFTEE} | tar -xf - -C $LOFTEE_DIR/
    ls $LOFTEE_DIR

    vep \\
    --vcf \\
    --format vcf \\
    -i {b.read_input(vcf_path)} \\
    --everything \\
    --allele_number \\
    --no_stats \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache $CACHE_DIR/vep/ \\
    --dir_plugins $LOFTEE_PLUGIN_PATH \\
    --fasta $FASTA \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())} \\
    -o {j.out_vcf} \\
    --compress_output bgzip
    """
    j.command(wrap_command(
        cmd, 
        setup_gcp=True, 
        monitor_space=True, 
    ))
    if out_vcf_path:
        b.write_output(j.out_vcf, out_vcf_path)
    return j
