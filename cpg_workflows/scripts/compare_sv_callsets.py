"""
test script, see what happens when we compare two of our own callsets
"""

from dataclasses import dataclass

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_references
)


@dataclass
class Dataset:
    name: str


truth_vcf = 'gs://cpg-seqr-test/gatk_sv/d84fe17c11fba10df5749878e62065f1624b35_3461/SpiceUpSVIDs/fresh_ids.vcf.bgz'
eval_vcf = 'gs://cpg-seqr-test/gatk_sv/b00828c6cf20df7160f20b23537148dbd8b270_2700/FilterGenotypes/filtered.vcf.gz'
input_dict: dict = {
    'output_prefix': 'concordance_test',
    'reference_dict': str(get_fasta().with_suffix('.dict')),
    'eval_vcf': eval_vcf,
    'truth_vcf': truth_vcf,
}
input_dict |= get_images(['gatk_docker', 'sv_base_mini_docker'])
input_dict |= get_references([{'contig_list': 'primary_contigs_list'}])

test_prefix = to_path('gs://cpg-seqr-test/gatk_sv/test_concordance')
expected_outs = {
    'concordance_vcf': test_prefix / 'sv_concordance.vcf.gz',
    'concordance_vcf_index': test_prefix / 'sv_concordance.vcf.gz.tbi',
}

# # can't use this in a rogue script... sigh
jobs = add_gatk_sv_jobs(
    dataset=Dataset('test_concordance'),
    wfl_name='SVConcordance',
    input_dict=input_dict,
    expected_out_dict=expected_outs
)
get_batch().run(wait=False)
