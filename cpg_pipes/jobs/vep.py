#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import images, utils, Path
from cpg_pipes.hb.batch import Batch
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.jobs import split_intervals
from cpg_pipes.jobs.vcf import gather_vcfs, subset_vcf
from cpg_pipes.refdata import RefData
from cpg_pipes.types import SequencingType


def vep(
    b: Batch,
    vcf_path: Path,
    refs: RefData,
    sequencing_type: SequencingType,
    out_vcf_path: Path | None = None,
    overwrite: bool = False,
    scatter_count: int | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF.
    """
    if out_vcf_path and utils.can_reuse(out_vcf_path, overwrite):
        return [b.new_job('VEP [reuse]', job_attrs)]

    scatter_count = scatter_count or RefData.number_of_joint_calling_intervals
    jobs: list[Job] = []
    intervals_j = split_intervals.get_intervals(
        b=b,
        refs=refs,
        sequencing_type=sequencing_type,
        scatter_count=scatter_count,
    )
    jobs.append(intervals_j)
    
    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    vep_jobs = []
    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        subset_j = subset_vcf(
            b,
            vcf=vcf,
            intervals=intervals_j[f'{idx}.interval_list'],
            refs=refs,
            job_attrs=(job_attrs or {}) | dict(intervals=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        j = _vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            refs=refs,
            job_attrs=(job_attrs or {}) | dict(intervals=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(j)
        vep_jobs.append(j)
    gather_j, gather_vcf = gather_vcfs(
        b=b,
        input_vcfs=[j.output_vcf['vcf.gz'] for j in vep_jobs],
        out_vcf_path=out_vcf_path,
    )
    jobs.append(gather_j)
    return jobs


def _vep_one(
    b: Batch,
    vcf: Path | hb.Resource,
    refs: RefData,
    out_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run VEP in a single job.
    """
    j = b.new_job('VEP', job_attrs)
    j.image(images.VEP_IMAGE)
    STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': '$LOFTEE_DIR/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': '$LOFTEE_DIR/human_ancestor.fa.gz',
        'conservation_file': '$LOFTEE_DIR/loftee.sql',
    }
    
    if not isinstance(vcf, hb.Resource):
        vcf = b.read_input(str(vcf))

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""\
    CACHE_DIR=/io/batch/cache
    LOFTEE_DIR=/io/batch/loftee
    LOFTEE_PLUGIN_PATH=/root/micromamba/share/ensembl-vep-105.0-0
    FASTA=$CACHE_DIR/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    mkdir -p $CACHE_DIR
    mkdir -p $LOFTEE_DIR
    gsutil cat {refs.vep_cache} | tar -xf - -C $CACHE_DIR/
    gsutil cat {refs.vep_loftee} | tar -xf - -C $LOFTEE_DIR/
    ls $LOFTEE_DIR

    vep \\
    --vcf \\
    --format vcf \\
    --compress_output bgzip \\
    -o {j.output_vcf['vcf.gz']} \\
    -i {vcf} \\
    --everything \\
    --allele_number \\
    --no_stats \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache $CACHE_DIR/vep/ \\
    --dir_plugins $LOFTEE_PLUGIN_PATH \\
    --fasta $FASTA \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())}
    
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """
    j.command(wrap_command(
        cmd, 
        setup_gcp=True, 
        monitor_space=True, 
        define_retry_function=True
    ))
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))
    return j
