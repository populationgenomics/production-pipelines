import logging
import os

import hail as hl

from cpg_utils import config, hail_batch, to_path
from talos import utils


COMPOSE_COMMAND = 'gcloud storage objects compose'

# one batch to rule them all
batch_instance = hail_batch.get_batch()


def compose_object_fragments(obj_folder: str, temp_dir: str, output_final: str):
    """
    Generate a gcloud compose command string to concatenate a list of VCF fragments
    Args:
        obj_folder (str): the GCS folder containing the VCF fragments

    Returns:
        the compose instructions to generate the output from the fragments
    """

    vcf_fragments = []
    with to_path(f'{obj_folder}/shard-manifest.txt').open() as read_handle:
        vcf_fragments = [os.path.join(obj_folder, line.strip()) for line in read_handle if line.strip()]

    if not vcf_fragments:
        raise ValueError('Input manifest is empty or contains no valid VCF fragments.')

    merge_round = 1
    condense_strings = []
    while len(vcf_fragments) > 1:
        new_fragments = []

        # If we have more fragments than the chunk size, we need to merge them in chunks
        for merge_chunk, fragment_list in enumerate(utils.chunks(vcf_fragments, 32)):
            output = f'{temp_dir}/{merge_round}/temp_chunk_{merge_chunk}.vcf.gz'
            condense_strings.append(f'{COMPOSE_COMMAND} {" ".join(fragment_list)} {output}')
            new_fragments.append(output)

        vcf_fragments = new_fragments
        merge_round += 1

    # one left, _move_ it to the non-tmp bucket
    # inter-bucket condense operations aren't valid, so we can't 'compose' from tmp to main
    condense_strings.append(f'gcloud storage mv {vcf_fragments[0]} {output_final}')

    # Write the final script to the output file
    script_path = os.path.join(temp_dir, 'compose_script.sh')
    with to_path(script_path).open('w') as script_file:
        script_file.write('#!/bin/bash\n\n')
        for condense_string in condense_strings:
            script_file.write(f'{condense_string}\n')
    return script_path


def export_and_condense(mt: hl.MatrixTable, sample: str, tmp_dir: str, output_final: str):
    """
    export a given MT as a VCF in parallel, and condense it into a single file
    """
    export_dir = os.path.join(tmp_dir, f'{sample}.vcf.bgz')

    shard_manifest = os.path.join(export_dir, 'shard-manifest.txt')

    if not to_path(shard_manifest).exists():
        logging.info(f'Creating {shard_manifest}')
        hl.export_vcf(mt, export_dir, parallel='separate_header')

    # create the composition commands
    command_script = compose_object_fragments(
        obj_folder=export_dir,
        temp_dir=os.path.join(tmp_dir, sample),
        output_final=output_final,
    )

    new_job = batch_instance.new_bash_job(f'Compose VCF for {sample}')
    new_job.image(config.config_retrieve(['workflow', 'driver_image']))

    localised_script = batch_instance.read_input(command_script)
    new_job.command(f'bash {localised_script}')


input_mt = 'gs://cpg-ghfm-kidgen-test/mt/bbf37e78fe_5866-ghfm-kidgen_full_copy.mt'
tmp_dir = 'gs://cpg-ghfm-kidgen-test-tmp/talos_benchmarking/ghfm-kidgen_ms_vcfs'
output_ms_vcfs = 'gs://cpg-ghfm-kidgen-test/talos_benchmarking/ghfm-kidgen_ms_vcfs'

logging.basicConfig(level=logging.INFO)

hail_batch.init_batch()

mt = hl.read_matrix_table(input_mt)

mt = mt.select_rows()

samples = sorted(mt.s.collect())

for chunk_num, sg_chunk in enumerate(utils.chunks(samples, 52)):
    out_path = os.path.join(output_ms_vcfs, f'chunk_{chunk_num}.vcf.bgz')
    sam_mt = mt.filter_cols(hl.literal(sg_chunk).contains(mt.s))
    sam_mt = sam_mt.filter_rows(hl.agg.any(sam_mt.GT.is_non_ref()))
    export_and_condense(sam_mt, f'chunk_{chunk_num}', tmp_dir, out_path)

batch_instance.run(wait=False)
