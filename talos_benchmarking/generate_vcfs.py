"""
Use the gs://cpg-acute-care-test/mt/bbf37e78fe_5866-acute-care_full_copy.mt

Using the whole range of variants, this script will grab all the samples in the MT and do the following:

- generate a subset of 250 samples, export their data as a multisample VCF, and single-sample VCFs
- of that 250, generate a subset of 100, export as a multisample VCF
- repeat for 50, 25, and 10 samples

These multi/single sample data will be the basis of Nextflow benchmarking for the Talos annotation pre-process
"""

import logging
import os
import random

import hail as hl

from cpg_utils import config, hail_batch, to_path
from cpg_workflows import utils

COMPOSE_COMMAND = 'gcloud storage objects compose'


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

    new_job = hail_batch.get_batch().new_bash_job(f'Compose VCF for {sample}')
    new_job.image(config.config_retrieve(['workflow', 'driver_image']))

    localised_script = hail_batch.get_batch().read_input(command_script)
    new_job.command(f'bash {localised_script}')


logging.basicConfig(level=logging.INFO)

random.seed(42)

# kick up a batch
hail_batch.init_batch()

# this was the original MT prior to generation of a smaller subset checkpoint
# input_mt = 'gs://cpg-acute-care-test/mt/bbf37e78fe_5866-acute-care_full_copy.mt'
# _250_samples = random.sample(samples, 250)

# this now contains only 250 samples, and limited to variants within that sample subset
input_mt = 'gs://cpg-acute-care-test/talos_benchmarking/250samples.mt'

output_prefix = 'gs://cpg-acute-care-test/talos_benchmarking/ms_vcfs'
ss_vcf_prefix = 'gs://cpg-acute-care-test/talos_benchmarking/solo_vcfs'

mt = hl.read_matrix_table(input_mt)

_250_samples = list(mt.s.collect())

tmp_dir = os.path.join(config.config_retrieve(['storage', 'default', 'tmp']), 'talos_benchmarking')


"""
replace the vcf extraction with https://github.com/populationgenomics/production-pipelines/pull/1270/files
"""

# for each_sam in _250_samples:
#     out_path = os.path.join(ss_vcf_prefix, f'{each_sam}.vcf.bgz')
#     if utils.exists(out_path):
#         logging.info(f'{out_path} exists, skipping')
#         continue
#     sam_mt = mt.filter_cols(mt.s == each_sam)
#     export_and_condense(sam_mt, each_sam, tmp_dir, out_path)

# export the bigboi
out_path = os.path.join(output_prefix, '250.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_250_samples).contains(mt.s))
    export_and_condense(sam_mt, '250_samples', tmp_dir, out_path)

# step down to 100 samples
_100_samples = random.sample(_250_samples, 100)
out_path = os.path.join(output_prefix, '100.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_100_samples).contains(mt.s))
    export_and_condense(sam_mt, '100_samples', tmp_dir, out_path)

# step down to 50 samples
_50_samples = random.sample(_100_samples, 50)
out_path = os.path.join(output_prefix, '50.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_50_samples).contains(mt.s))
    export_and_condense(sam_mt, '50_samples', tmp_dir, out_path)

# step down to 25 samples
_25_samples = random.sample(_50_samples, 25)
out_path = os.path.join(output_prefix, '25.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_25_samples).contains(mt.s))
    export_and_condense(sam_mt, '25_samples', tmp_dir, out_path)

# step down to 10 samples
_10_samples = random.sample(_25_samples, 10)
out_path = os.path.join(output_prefix, '10.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_10_samples).contains(mt.s))
    export_and_condense(sam_mt, '10_samples', tmp_dir, out_path)

# step down to 5 samples
_5_samples = random.sample(_10_samples, 5)
out_path = os.path.join(output_prefix, '5.vcf.bgz')
if not utils.exists(out_path):
    logging.info(f'Creating {out_path}')
    sam_mt = mt.filter_cols(hl.literal(_5_samples).contains(mt.s))
    export_and_condense(sam_mt, '5_samples', tmp_dir, out_path)

hail_batch.get_batch().run(wait=False)