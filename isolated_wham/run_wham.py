#!/usr/bin/env python3

"""
Runs Wham on a CRAM file
"""

from typing import Any

from cpg_utils import to_path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    CromwellJobSizes,
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references,
)


class WhamDataset:
    name = 'wham_test'


cram = 'gs://cpg-thousand-genomes-test/cram/CPG196519.cram'

input_dict: dict[str, Any] = dict(
    bam_or_cram_file=cram,
    bam_or_cram_index=f'{cram}.crai',
    sample_id='HG00096',
    reference_fasta=str(get_fasta()),
    reference_index=str(get_fasta()) + '.fai',
)
input_dict |= get_images(['wham_docker'])
input_dict |= get_references(
    [
        'primary_contigs_list',
        {'include_bed_file': 'wham_include_list_bed_file'},
    ],
)

expected_out = {
    'vcf': to_path(output_path('CPG196519.vcf.gz')),
    'index': to_path(output_path('CPG196519.vcf.gz.tbi'))
}


add_gatk_sv_jobs(
    dataset=WhamDataset(),
    wfl_name='Whamg',
    input_dict=input_dict,
    expected_out_dict=expected_out,
    sequencing_group_id='HG00096',
    job_size=CromwellJobSizes.LARGE,
)

get_batch().run(wait=False)
