"""
concept test of a stage using a test_mode
"""
import inspect
import re
from cpg_workflows.stages.happy_validation import ValidationMtToVcf, ValidationHappyOnVcf
from cpg_workflows.query_modules.validation import single_sample_vcf_from_dataset_vcf
from unittest import mock

from cpg_workflows.batch import get_batch

from .. import set_config

from test.test_workflow import mock_create_cohort


TOML = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = 'test'
sequencing_type = 'genome'
only_stages = ["{only}"]

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

[storage.default]
default = "{directory}"

[storage.fewgenomes]
default = "{directory}"

[storage.my_dataset]
default = "{directory}"

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'

[references]
refgenome_sdf = "refgenome_sdf"
stratification = "{directory}"

[references.broad]
ref_fasta = "ref_fasta"

[references.SAMPLE1]
vcf = "vcf"
index = "index"
bed = "bed"

[references.SAMPLE2]
vcf = "vcf"
index = "index"
bed = "bed"

[inputs]
sample_hash = "boop"

[images]
cpg_workflows = "test"
hap-py = "test"
"""


@mock.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)
@mock.patch('cpg_workflows.workflow.list_of_all_dir_contents', set)
def test_validation_mt_to_vcf(tmp_path):

    from cpg_workflows.workflow import run_workflow

    conf = TOML.format(directory=tmp_path, only='ValidationMtToVcf')
    set_config(conf, tmp_path / 'config.toml')
    run_workflow(stages=[ValidationMtToVcf], test_mode=True)
    jobs = get_batch()._jobs

    expected_code = inspect.getsource(single_sample_vcf_from_dataset_vcf)
    assert len(jobs) == 2
    for job in jobs:
        assert job.attributes['stage'] == 'ValidationMtToVcf'
        assert expected_code in job._command[0]


@mock.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)
@mock.patch('cpg_workflows.workflow.list_of_all_dir_contents', set)
def test_validation_happy_on_vcf(tmp_path):

    from cpg_workflows.workflow import run_workflow

    conf = TOML.format(directory=tmp_path, only='ValidationHappyOnVcf')
    # create the definition file we check for
    definition_path = tmp_path / 'definition.tsv'
    definition_path.write_text('test')

    set_config(conf, tmp_path / 'config.toml')
    run_workflow(stages=[ValidationHappyOnVcf], test_mode=True)
    jobs = get_batch()._jobs

    assert len(jobs) == 2
    for job in jobs:
        j_dict = job.__dict__
        token = j_dict['_token']
        sg = j_dict['attributes']['sequencing_group']
        participant = j_dict['attributes']['participant_id']
        whole_bit = re.compile(
            f'mkdir \${{BATCH_TMPDIR}}/Run_Happy_on_{participant}_VCF-{token}/output && '
            f'hap.py \${{BATCH_TMPDIR}}/inputs/\w+/vcf \$\{{BATCH_TMPDIR}}/inputs/\w+/{sg}.vcf.bgz '
            '-r \${BATCH_TMPDIR}/inputs/\w+/ref_fasta '
            '-R \${BATCH_TMPDIR}/inputs/\w+/bed '
            f'-o \${{BATCH_TMPDIR}}/Run_Happy_on_{participant}_VCF-{token}/output/output '
            '--leftshift '
            '--threads 4 '
            '--preprocess-truth '
            '--engine-vcfeval-path=/opt/hap.py/libexec/rtg-tools-install/rtg '
            '--engine-vcfeval-template \${BATCH_TMPDIR}/inputs/\w+ '
            '--engine=vcfeval '
            '--stratification \${BATCH_TMPDIR}/inputs/\w+/definition.tsv'
        )
        assert job.attributes['stage'] == 'ValidationHappyOnVcf'
        assert re.match(whole_bit, job._command[0])
