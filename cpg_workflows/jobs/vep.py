#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""
import logging
from typing import Literal

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.jobs.vcf import gather_vcfs
from cpg_workflows.query_modules import vep
from cpg_workflows.utils import can_reuse


def add_vep_jobs(
    tmp_prefix: Path,
    input_vcf_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.

    Changed structure - this no longer does interval splitting on an input file, if you want a divided VCF, do that
    elsewhere. The pipeline uses this in one way - pass a list of sub-VCFs, annotate each, re-join.
    In production there's no need for this to also know how to split
    """
    to_hail_table = out_path and out_path.suffix == '.ht'

    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and can_reuse(out_path):
        return []

    jobs: list[Job] = []

    input_vcf_parts: list[hb.ResourceGroup] = [
        get_batch().read_input_group(**{'vcf.gz': str(path), 'vcf.gz.tbi': str(path) + '.tbi'})
        for path in input_vcf_paths
    ]

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []
    for idx, resource in enumerate(input_vcf_parts, 1):
        if to_hail_table:
            result_part_path = result_parts_bucket / f'part{idx}.jsonl'
        else:
            result_part_path = result_parts_bucket / f'part{idx}.vcf.gz'
        result_part_paths.append(result_part_path)
        if can_reuse(result_part_path):
            continue

        vep_one_job = vep_one(
            vcf=resource['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=result_part_paths[idx],
            job_attrs=(job_attrs or {}) | dict(part=f'{idx}'),
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    if to_hail_table:
        j = gather_vep_json_to_ht(
            vep_results_paths=result_part_paths,
            out_path=out_path,
            job_attrs=job_attrs,
        )
        if jobs:
            j.depends_on(*jobs)
        jobs.append(j)

    elif len(result_part_paths) != 1:
        gather_jobs = gather_vcfs(
            b=get_batch(),
            input_vcfs=result_part_paths,
            out_vcf_path=out_path,
        )
        if gather_jobs and jobs:
            for j in gather_jobs:
                j.depends_on(*jobs)
        jobs.extend(gather_jobs)
    return jobs


def gather_vep_json_to_ht(
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """

    vep_version = get_config()['workflow']['vep_version']

    j = get_batch().new_job('VEP', job_attrs)
    j.image(image_path('cpg_workflows'))
    j.command(
        query_command(
            vep,
            vep.vep_json_to_ht.__name__,
            [str(p) for p in vep_results_paths],
            str(out_path),
            vep_version,
            setup_gcp=True,
        ),
    )
    return j


def vep_one(
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

    vep_version = get_config()['workflow'].get('vep_version')
    if not vep_version:
        raise IndexError('No VEP version specified in config.workflow.vep_version')

    # check that the cache and image for this version exist
    vep_image = image_path(f'vep_{vep_version}')
    vep_mount_path = to_path(reference_path(f'vep_{vep_version}_mount'))
    assert all([vep_image, vep_mount_path])
    logging.info(f'Using VEP {vep_version}')

    j = get_batch().new_job('VEP', (job_attrs or {}) | dict(tool=f'VEP {vep_version}'))
    j.image(vep_image)
    splice_ai = get_config()['workflow'].get('spliceai_plugin', False)

    # vep is single threaded, with a middling memory requirement
    # during test it can exceed 8GB, so we'll give it 16GB
    j.memory('16Gi').storage('15Gi').cpu(1)

    if not isinstance(vcf, hb.ResourceFile):
        vcf = get_batch().read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
        assert isinstance(j.output, hb.ResourceGroup)
        output = j.output['vcf.gz']
    else:
        assert isinstance(j.output, hb.ResourceFile)
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    # assume for now that only VEP 105 has a non-standard install location
    loftee_conf = {
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
        'loftee_path': ('$MAMBA_ROOT_PREFIX/share/ensembl-vep' if vep_version == '105' else '$VEP_DIR_PLUGINS'),
    }

    # sexy new plugin - only present in 110 build
    alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '

    # UTRannotator plugin doesn't support JSON output at this time; only activate for VCF outputs
    # VCF annotation doesn't utilise the aggregated Seqr reference data, including spliceAI
    # SpliceAI requires both indel and SNV files to be present (~100GB), untested
    vcf_plugins = '--plugin UTRAnnotator,file=$UTR38 '
    if splice_ai:
        vcf_plugins += (
            f'--plugin SpliceAI,snv={vep_dir}/spliceai_scores.raw.snv.hg38.vcf.gz,'
            f'indel={vep_dir}/spliceai_scores.raw.indel.hg38.vcf.gz '
        )

    # VEP 105 installs plugins in non-standard locations
    loftee_plugin_path = '--dir_plugins $MAMBA_ROOT_PREFIX/share/ensembl-vep '

    cmd = f"""\
    FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz
    vep \\
    --format vcf \\
    --{out_format} {'--compress_output bgzip' if out_format == 'vcf' else ''} \\
    -o {output} \\
    -i {vcf} \\
    --everything \\
    --mane_select \\
    --allele_number \\
    --minimal \\
    --species homo_sapiens \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --fasta $FASTA \\
    {loftee_plugin_path if vep_version == '105' else alpha_missense_plugin} \
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())} \
    {vcf_plugins if (vep_version == '110' and out_format == 'vcf') else ''}
    """

    if out_format == 'vcf':
        cmd += f'tabix -p vcf {output}'

    j.command(cmd)

    if out_path:
        get_batch().write_output(j.output, str(out_path).replace('.vcf.gz', ''))

    return j
