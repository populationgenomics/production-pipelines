#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""

from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import (
    image_path,
    reference_path,
    command,
    authenticate_cloud_credentials_in_job,
)
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse

from .picard import get_intervals
from .vcf import gather_vcfs, subset_vcf


def add_vep_jobs(
    b: Batch,
    vcf_path: Path,
    tmp_prefix: Path,
    scatter_count: int,
    out_path: Path | None = None,
    overwrite: bool = False,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and can_reuse(out_path, overwrite):
        return []

    jobs: list[Job] = []
    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    parts_bucket = tmp_prefix / 'vep' / 'parts'
    part_files = []

    intervals_j, intervals = get_intervals(
        b=b,
        scatter_count=scatter_count,
        job_attrs=job_attrs,
        output_prefix=tmp_prefix / f'intervals_{scatter_count}',
    )
    if intervals_j:
        jobs.append(intervals_j)

    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        if to_hail_table:
            part_path = parts_bucket / f'part{idx + 1}.jsonl'
        else:
            part_path = parts_bucket / f'part{idx + 1}.vcf.gz'
        part_files.append(part_path)
        if can_reuse(part_path):
            continue

        subset_j = subset_vcf(
            b,
            vcf=vcf,
            interval=intervals[idx],
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        assert isinstance(subset_j.output_vcf, hb.ResourceGroup)

        # noinspection PyTypeChecker
        vep_one_job = vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=part_path,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            overwrite=overwrite,
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    gather_j: Job | None = None
    if to_hail_table:
        gather_j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=part_files,
            out_path=out_path,
            job_attrs=job_attrs,
            depends_on=jobs,
        )
    else:
        assert len(part_files) == scatter_count
        gather_j, gather_vcf = gather_vcfs(
            b=b,
            input_vcfs=part_files,
            out_vcf_path=out_path,
        )
    if gather_j:
        gather_j.depends_on(*jobs)
        jobs.append(gather_j)
    return jobs


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
    depends_on: list[hb.job.Job] | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    # Importing this requires CPG_CONFIG_PATH to be already set, that's why
    # we are not importing it on the top level.
    from analysis_runner import dataproc

    script_path = to_path(__file__).parent / 'dataproc_scripts' / 'vep_json_to_ht.py'

    j = dataproc.hail_dataproc_job(
        b,
        f'{script_path} --out-path {out_path} '
        + ' '.join(str(p) for p in vep_results_paths),
        max_age='24h',
        packages=[
            'cpg_workflows',
            'google',
            'fsspec',
            'gcloud',
        ],
        num_workers=2,
        num_secondary_workers=20,
        job_name=f'VEP JSON to Hail table',
        depends_on=depends_on,
        scopes=['cloud-platform'],
        pyfiles=['query_modules'],
        init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
    )
    j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
    return j


def vep_one(
    b: Batch,
    vcf: Path | hb.ResourceFile,
    out_path: Path | None = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run a single VEP job.
    """
    if out_path and can_reuse(out_path, overwrite):
        return None

    j = b.new_job('VEP', (job_attrs or {}) | dict(tool='vep'))
    j.image(image_path('vep'))
    STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    if not isinstance(vcf, hb.ResourceFile):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(
            output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )
        assert isinstance(j.output, hb.ResourceGroup)
        output = j.output['vcf.gz']
    else:
        assert isinstance(j.output, hb.ResourceFile)
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    vep_mount_path = reference_path('vep_mount')
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])
    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
    }

    authenticate_cloud_credentials_in_job(j)
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
        command(
            cmd,
            setup_gcp=True,
            monitor_space=True,
        )
    )
    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))
    return j
