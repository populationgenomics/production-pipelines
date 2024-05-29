#! /usr/bin/env python3

"""
this demo script takes a long-read SV VCF, defined in config, and annotates it
"""
from cpg_utils import to_path
from cpg_utils.config import config_retrieve, output_path, try_get_ar_guid
from cpg_utils.hail_batch import get_batch

from cpg_workflows.stages.gatk_sv.gatk_sv_common import get_references, get_images, add_gatk_sv_jobs


class DataDude:
    name: str = 'rdnow'


if __name__ == '__main__':
    output_folder = to_path(output_path('sv_test_annotate'))
    expected_out = {
        'annotated_vcf': output_folder / 'filtered_annotated.vcf.bgz',
        'annotated_vcf_index': output_folder / 'filtered_annotated.vcf.bgz.tbi',
    }

    input_vcf = to_path(config_retrieve(['workflow', 'input_vcf']))

    input_dict: dict = {
        'vcf': input_vcf,
        'prefix': 'rdnow',
        'ped_file': to_path(config_retrieve(['workflow', 'ped_file'])),
        'sv_per_shard': 5000,
        'population': config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'ref_prefix': config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'use_hail': False,
    }

    input_dict |= get_references(
        [
            'noncoding_bed',
            'protein_coding_gtf',
            {'ref_bed': 'external_af_ref_bed'},
            {'contig_list': 'primary_contigs_list'},
        ],
    )

    # images!
    input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker', 'gatk_docker'])

    jobs = add_gatk_sv_jobs(
        dataset=DataDude(),
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=expected_out,
        labels={'ar-guid': try_get_ar_guid()},
    )

    get_batch().run(wait=False)
