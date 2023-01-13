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
    query_command,
)
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse

from .picard import get_intervals
from .vcf import gather_vcfs, subset_vcf


def add_vep_jobs(
    b: Batch,
    input_siteonly_vcf_path: Path,
    tmp_prefix: Path,
    scatter_count: int,
    input_siteonly_vcf_part_paths: list[Path] | None = None,
    out_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and can_reuse(out_path):
        return []

    jobs: list[Job] = []
    siteonly_vcf = b.read_input_group(
        **{
            'vcf.gz': str(input_siteonly_vcf_path),
            'vcf.gz.tbi': str(input_siteonly_vcf_path) + '.tbi',
        }
    )

    input_vcf_parts: list[hb.ResourceGroup] = []
    if input_siteonly_vcf_part_paths:
        assert len(input_siteonly_vcf_part_paths) == scatter_count
        for path in input_siteonly_vcf_part_paths:
            input_vcf_parts.append(
                b.read_input_group(
                    **{'vcf.gz': str(path), 'vcf.gz.tbi': str(path) + '.tbi'}
                )
            )
    else:
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
            subset_j = subset_vcf(
                b,
                vcf=siteonly_vcf,
                interval=intervals[idx],
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            jobs.append(subset_j)
            assert isinstance(subset_j.output_vcf, hb.ResourceGroup)
            input_vcf_parts.append(subset_j.output_vcf)

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []
    for idx in range(scatter_count):
        if to_hail_table:
            result_part_path = result_parts_bucket / f'part{idx + 1}.jsonl'
        else:
            result_part_path = result_parts_bucket / f'part{idx + 1}.vcf.gz'
        result_part_paths.append(result_part_path)
        if can_reuse(result_part_path):
            continue

        # noinspection PyTypeChecker
        vep_one_job = vep_one(
            b,
            vcf=input_vcf_parts[idx]['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=result_part_paths[idx],
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    if to_hail_table:
        j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=result_part_paths,
            out_path=out_path,
            job_attrs=job_attrs,
            depends_on=jobs,
        )
        gather_jobs = [j]
    else:
        assert len(result_part_paths) == scatter_count
        gather_jobs = gather_vcfs(
            b=b,
            input_vcfs=result_part_paths,
            out_vcf_path=out_path,
        )
    for j in gather_jobs:
        j.depends_on(*jobs)
        jobs.append(j)
    return jobs


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
    use_dataproc: bool = False,
    depends_on: list[hb.job.Job] | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    if use_dataproc:
        # Importing this requires CPG_CONFIG_PATH to be already set, that's why
        # we are not importing it on the top level.
        from analysis_runner import dataproc

        # Script path and pyfiles should be relative to the repository root
        script_path = 'cpg_workflows/dataproc_scripts/vep_json_to_ht.py'
        pyfiles = ['cpg_workflows/query_modules']

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
            pyfiles=pyfiles,
            init=['gs://cpg-common-main/hail_dataproc/install_common.sh'],
        )
        j.attributes = (job_attrs or {}) | {'tool': 'hailctl dataproc'}
    else:
        from cpg_workflows.query_modules import vep

        j = b.new_job(f'VEP', job_attrs)
        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                vep,
                vep.vep_json_to_ht.__name__,
                [str(p) for p in vep_results_paths],
                str(out_path),
                setup_gcp=True,
            )
        )
    return j


def vep_one(
    b: Batch,
    vcf: Path | hb.ResourceFile,
    out_path: Path | None = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: dict | None = None,
) -> Job | None:
    """
    Run a single VEP job.
    """
    if out_path and can_reuse(out_path):
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

    LOFTEE_PLUGIN_PATH=$MAMBA_ROOT_PREFIX/share/ensembl-vep
    FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz

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
