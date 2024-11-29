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
from cpg_workflows.utils import generator_chunks, check_exists_path

VEP_CHUNK_SIZE = config_retrieve(['AnnotateFragmentedVcfWithVep', 'vep_chunk_size'], 30)


def add_vep_jobs(
    input_vcfs: list[ResourceGroup],
    tmp_prefix: Path,
    final_out_path: Path,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes annotations as JSON output

    More experimental batching

    Args:
        input_vcfs ():
        tmp_prefix ():
        final_out_path ():
        job_attrs ():
    """
    # check that the cache and image for this version exist
    vep_image = image_path('vep_110')
    vep_mount_path = to_path(reference_path('vep_110_mount'))
    assert all([vep_image, vep_mount_path])

    jobs: list[Job] = []

    fragment_count = len(input_vcfs)

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'

    result_part_paths = [result_parts_bucket / f'part{idx + 1}.jsonl' for idx in range(fragment_count)]
    for chunk_number, chunk_data in enumerate(generator_chunks(zip(input_vcfs, result_part_paths), VEP_CHUNK_SIZE), 1):

        j = get_batch().new_job(f'VEP part {chunk_number}', (job_attrs or {}) | dict(tool='VEP 110'))
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

        alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '

        for idx, (vcf, out_path) in enumerate(chunk_data):
            if check_exists_path(out_path):
                continue

            # string index value, used to separate the various outputs
            stringdex = str(idx)

            j.command(f"""
            set -x
            vep \\
            --format vcf \\
            --json \\
            -o {j[stringdex]} \\
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
            """)

            get_batch().write_output(j[stringdex], str(out_path))

    gather_job = get_batch().new_job('VEP Gather', job_attrs)
    gather_job.depends_on(*jobs)

    gather_job.image(config_retrieve(['workflow', 'driver_image']))
    json_paths = ' '.join(str(p) for p in result_part_paths)
    gather_job.command(f'vep_json_to_ht --input {json_paths} --output {final_out_path}')
    jobs.append(gather_job)
    return jobs
