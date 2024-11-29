#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""

from textwrap import dedent

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch.resource import ResourceGroup

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.utils import VCF_GZ, can_reuse


def add_vep_jobs(
    input_vcfs: list[ResourceGroup],
    tmp_prefix: Path,
    final_out_path: Path,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes annotations as JSON output

    Args:
        input_vcfs (list[ResourceGroup]): List of input VCFs, localised into the Hail Batch
        tmp_prefix (Path): Path to the temporary directory for writing fragments of annotation output
        final_out_path (Path): Path to write the final annotation output
        job_attrs (dict | None): Job attributes for the Hail Batch job
    """

    jobs: list[Job] = []

    fragment_count = len(input_vcfs)

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []
    for idx, resource in enumerate(input_vcfs):
        result_part_path = result_parts_bucket / f'part{idx + 1}.jsonl'

        result_part_paths.append(result_part_path)

        if can_reuse(result_part_path):
            continue

        jobs.append(
            vep_one(
                vcf=resource[VCF_GZ],
                out_path=str(result_part_paths[idx]),
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{fragment_count}'),
            ),
        )

    j = gather_vep_json_to_ht(
        vep_results_paths=result_part_paths,
        out_path=final_out_path,
        job_attrs=job_attrs,
    )
    j.depends_on(*jobs)
    jobs.append(j)
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

    j = get_batch().new_job('VEP', job_attrs)
    j.image(config_retrieve(['workflow', 'driver_image']))
    json_paths = ' '.join(str(p) for p in vep_results_paths)
    j.command(f'vep_json_to_ht --input {json_paths} --output {out_path}')
    return j


def vep_one(
    vcf: Path | hb.ResourceFile,
    out_path: str,
    job_attrs: dict | None = None,
) -> Job:
    """
    Run a single VEP job.

    Args:
        vcf ():
        out_path ():
        job_attrs ():
    """

    # check that the cache and image for this version exist
    vep_image = image_path('vep_110')
    vep_mount_path = to_path(reference_path('vep_110_mount'))
    assert all([vep_image, vep_mount_path])

    j = get_batch().new_job('VEP', (job_attrs or {}) | dict(tool='VEP 110'))
    j.image(vep_image)

    # vep is single threaded, with a middling memory requirement
    # during test it can exceed 8GB, so we'll give it 16GB
    j.memory('16Gi').storage('15Gi').cpu(1)

    # gcsfuse works only with the root bucket, without prefix:
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    # assume for now that only VEP 105 has a non-standard install location
    loftee_conf = {
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
        'loftee_path': '$VEP_DIR_PLUGINS',
    }

    # sexy new plugin - only present in 110 build
    alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '

    j.command(
        dedent(
            f"""
    set -x
    vep \\
    --format vcf \\
    --json \\
    -o {j.output} \\
    -i {vcf} \\
    --everything \\
    --mane_select \\
    --allele_number \\
    --minimal \\
    --species homo_sapiens \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --fasta /cpg-common-main/references/vep/110/mount/vep/homo_sapiens/110/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
    {alpha_missense_plugin} \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())} \\
    --plugin UTRAnnotator,file=$UTR38
    """,
        ),
    )

    get_batch().write_output(j.output, out_path)

    return j
