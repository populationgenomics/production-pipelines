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


cram = 'gs://cpg-seqr-test/thousand_genomes_copy/HG00096.cram'

input_dict: dict[str, Any] = dict(
    bam_or_cram_file=cram,
    bam_or_cram_index=f'{cram}.crai',
    sample_id='HG00096',
    reference_fasta=str(get_fasta()),
    reference_index=str(get_fasta()) + '.fai',
    reference_dict=str(get_fasta().with_suffix('.dict')),
    reference_version='38',
)
input_dict |= get_images(['wham_docker'])
input_dict |= get_references(
    [
        'primary_contigs_list',
        {'include_bed_file': 'wham_include_list_bed_file'},
        {'sd_locs_vcf': 'dbsnp_vcf'},
    ],
)

expected_out = {
    'vcf': to_path(output_path('HG00096.vcf.gz')),
    'index': to_path(output_path('HG00096.vcf.gz.tbi'))
}


add_gatk_sv_jobs(
    dataset=WhamDataset(),
    wfl_name='Whamg.wdl',
    input_dict=input_dict,
    expected_out_dict=expected_out,
    sequencing_group_id='HG00096',
    job_size=CromwellJobSizes.LARGE,
)

get_batch().run(wait=False)
