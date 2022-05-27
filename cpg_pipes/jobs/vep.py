#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""
import logging
from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import utils, Path, to_path
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.hb.command import wrap_command, python_command
from cpg_pipes.jobs import split_intervals
from cpg_pipes.jobs.vcf import gather_vcfs, subset_vcf
from cpg_pipes.providers.images import Images
from cpg_pipes.query import vep as vep_module
from cpg_pipes.providers.refdata import RefData
from cpg_pipes.types import SequencingType


logger = logging.getLogger(__file__)


def vep_jobs(
    b: Batch,
    vcf_path: Path,
    refs: RefData,
    images: Images,
    hail_billing_project: str,
    hail_bucket: Path,
    tmp_bucket: Path,
    out_path: Path | None = None,
    overwrite: bool = False,
    scatter_count: int | None = RefData.number_of_vep_intervals,
    sequencing_type: SequencingType = SequencingType.GENOME,
    intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Whites a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and str(out_path).rstrip('/').endswith('.ht')
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and utils.can_reuse(out_path, overwrite):
        return [b.new_job('VEP [reuse]', job_attrs)]

    scatter_count = scatter_count or RefData.number_of_vep_intervals
    jobs: list[Job] = []
    intervals_j, intervals = split_intervals.get_intervals(
        b=b,
        refs=refs,
        images=images,
        sequencing_type=sequencing_type,
        intervals_path=intervals_path,
        scatter_count=scatter_count,
    )
    jobs.append(intervals_j)

    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    parts_bucket = tmp_bucket / 'vep' / 'parts'
    part_files = []

    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        subset_j = subset_vcf(
            b,
            vcf=vcf,
            interval=intervals[idx],
            refs=refs,
            images=images,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        if to_hail_table:
            part_path = parts_bucket / f'part{idx + 1}.json_list'
        else:
            part_path = None
        # noinspection PyTypeChecker
        j = vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=part_path,
            refs=refs,
            images=images,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            overwrite=overwrite,
        )
        jobs.append(j)
        if to_hail_table:
            part_files.append(part_path)
        else:
            part_files.append(j.output['vcf.gz'])

    if to_hail_table:
        gather_j = gather_vep_json_to_ht(
            b=b,
            images=images,
            vep_results_paths=part_files,
            hail_billing_project=hail_billing_project,
            hail_bucket=hail_bucket,
            out_path=out_path,
            job_attrs=job_attrs,
        )
    else:
        assert len(part_files) == scatter_count
        gather_j, gather_vcf = gather_vcfs(
            b=b,
            images=images,
            input_vcfs=part_files,
            out_vcf_path=out_path,
        )
    gather_j.depends_on(*jobs)
    jobs.append(gather_j)
    return jobs


def gather_vep_json_to_ht(
    b: Batch,
    images: Images,
    vep_results_paths: list[Path],
    hail_billing_project: str,
    hail_bucket: Path,
    out_path: Path,
    job_attrs: dict | None = None,
):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    j = b.new_job('VEP json to Hail table', job_attrs)
    j.image(images.get('hail'))
    cmd = python_command(
        vep_module,
        vep_module.vep_json_to_ht.__name__,
        [str(p) for p in vep_results_paths],
        str(out_path),
        setup_gcp=True,
        hail_billing_project=hail_billing_project,
        hail_bucket=str(hail_bucket),
        default_reference=RefData.genome_build,
    )
    j.command(cmd)
    return j


def vep_one(
    b: Batch,
    vcf: Path | hb.Resource,
    refs: RefData,
    images: Images,
    out_path: Path | None = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Run a single VEP job.
    """
    j = b.new_job('VEP', job_attrs)
    if out_path and utils.can_reuse(out_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.get('vep'))
    STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    if not isinstance(vcf, hb.Resource):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(
            output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )
        output = j.output['vcf.gz']
    else:
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    base_bucket_name = refs.vep_mount.drive
    data_mount = to_path(f'/{base_bucket_name}')
    j.cloudfuse(base_bucket_name, str(data_mount), read_only=True)
    vep_dir = data_mount / 'vep' / 'GRCh38'
    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
    }

    cmd = f"""\
    ls {vep_dir}
    ls {vep_dir}/vep

    LOFTEE_PLUGIN_PATH=/root/micromamba/share/ensembl-vep-105.0-1
    FASTA={vep_dir}/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    vep \\
    --format vcf \\
    --{out_format} {'--compress_output bgzip' if out_format == 'vcf' else ''} \\
    -o {output} \\
    -i {vcf} \\
    --everything \\
    --allele_number \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --dir_plugins $LOFTEE_PLUGIN_PATH \\
    --fasta $FASTA \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())}
    """
    if out_format == 'vcf':
        cmd += f'tabix -p vcf {output}'

    j.command(
        wrap_command(
            cmd,
            setup_gcp=True,
            monitor_space=True,
        )
    )
    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))
    return j
