"""
Test seqr-loader workflow.
"""
from unittest.mock import mock_open

import toml
from cpg_utils import to_path
from pytest_mock import MockFixture
from . import results_prefix, update_dict

TOML = f"""
[workflow]
dataset_gcp_project = "test-analysis-dataset-1234"
dataset = "test-analysis-dataset"
access_level = "test"
sequencing_type = "genome"
driver_image = "<stub>"
skip_stages = [ "Align",]
check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = "local"
local_dir = "{results_prefix}"

[hail]
billing_project = "test-analysis-dataset"
delete_scratch_on_exit = true
dry_run = true
backend = "local"

[images]
bcftools = "stub"
bedtools = "stub"
fastqc = "stub"
gatk = "stub"
hap-py = "stub"
multipy = "stub"
multiqc = "stub"
peer = "stub"
picard = "stub"
samtools = "stub"
somalier = "stub"
vep = "stub"
verifybamid = "stub"

[references]
genome_build = "GRCh38"
vep_mount = "stub"
liftover_38_to_37 = "stub"
somalier_sites = "stub"
seqr_combined_reference_data = "stub"
seqr_clinvar = "stub"

[storage.default]
default = '{results_prefix()}'
web = "{results_prefix()}-web"
analysis = "{results_prefix()}-analysis"
tmp = "{results_prefix()}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.test-input-dataset]
default = "{results_prefix()}"
web = "{results_prefix()}-web"
analysis = "{results_prefix()}-analysis"
tmp = "{results_prefix()}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.test-analysis-dataset]
default = "{results_prefix()}"
web = "{results_prefix()}-web"
analysis = "{results_prefix()}-analysis"
tmp = "{results_prefix()}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[references.broad]
dragmap_prefix = "stub"
ref_fasta = "stub"
noalt_bed = "stub"
genome_calling_interval_lists = "stub"
exome_calling_interval_lists = "stub"
genome_evaluation_interval_lists = "stub"
exome_evaluation_interval_lists = "stub"
genome_coverage_interval_list = "stub"
unpadded_intervals_file = "stub"
dbsnp_vcf = "stub"
dbsnp_vcf_index = "stub"
hapmap_vcf = "stub"
hapmap_vcf_index = "stub"
omni_vcf = "stub"
omni_vcf_index = "stub"
one_thousand_genomes_vcf = "stub"
one_thousand_genomes_vcf_index = "stub"
mills_vcf = "stub"
mills_vcf_index = "stub"
axiom_poly_vcf = "stub"
axiom_poly_vcf_index = "stub"
genome_contam_ud = "stub"
genome_contam_bed = "stub"
genome_contam_mu = "stub"
exome_contam_ud = "stub"
exome_contam_bed = "stub"
exome_contam_mu = "stub"

[references.gnomad]
tel_and_cent_ht = "stub"
lcr_intervals_ht = "stub"
seg_dup_intervals_ht = "stub"
clinvar_ht = "stub"
hapmap_ht = "stub"
kgp_omni_ht = "stub"
kgp_hc_ht = "stub"
mills_ht = "stub"
"""


def _mock_config() -> dict:
    d: dict = {}
    for fp in [
        to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml',
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'seqr_loader.toml',
    ]:
        with fp.open():
            update_dict(d, toml.load(fp))

    update_dict(d, toml.loads(TOML))
    return d


def _mock_cohort():
    from cpg_workflows.targets import Cohort
    from cpg_workflows.filetypes import BamPath, FastqPair, FastqPairs

    cohort = Cohort()
    ds = cohort.create_dataset('test-input-dataset')
    ds.add_sample(
        'CPG01',
        'SAMPLE1',
        alignment_input_by_seq_type={
            'genome': BamPath('gs://test-input-dataset-upload/sample1.bam')
        },
    )
    ds.add_sample(
        'CPG02',
        'SAMPLE2',
        alignment_input_by_seq_type={
            'genome': FastqPairs(
                [
                    FastqPair(
                        'gs://test-input-dataset-upload/sample2_L1_R1.fq.gz',
                        'gs://test-input-dataset-upload/sample2_L1_R2.fq.gz',
                    ),
                    FastqPair(
                        'gs://test-input-dataset-upload/sample2_L2_R1.fq.gz',
                        'gs://test-input-dataset-upload/sample2_L2_R2.fq.gz',
                    ),
                ]
            )
        },
    )
    return cohort


def test_seqr_loader_dry(mocker: MockFixture):
    """
    Test entire seqr-loader in a dry mode.
    """
    mocker.patch('cpg_utils.config.get_config', _mock_config)
    mocker.patch('cpg_workflows.inputs.create_cohort', _mock_cohort)

    def mock_exists(*args, **kwargs) -> bool:
        return False

    def do_nothing(*args, **kwargs):
        return None

    def mock_create_new_analysis(*args, **kwargs) -> int:
        return 1

    mocker.patch('pathlib.Path.open', mock_open(read_data='<stub>'))
    # functions like get_intervals checks file existence
    mocker.patch('cloudpathlib.cloudpath.CloudPath.exists', mock_exists)
    # cloudfuse (used in Vep) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.cloudfuse', do_nothing)
    # always_run (used in MtToEs -> hail_dataproc_job) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.always_run', do_nothing)
    # can't access secrets from CI environment
    mocker.patch(
        'cpg_workflows.stages.seqr_loader.es_password', lambda: 'test-password'
    )
    mocker.patch(
        'sample_metadata.apis.AnalysisApi.create_new_analysis',
        mock_create_new_analysis,
    )
    mocker.patch('sample_metadata.apis.AnalysisApi.update_analysis_status', do_nothing)

    from cpg_workflows.batch import get_batch
    from cpg_workflows.inputs import get_cohort
    from cpg_workflows.stages.cram_qc import CramMultiQC
    from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
    from cpg_workflows.stages.joint_genotyping_qc import JointVcfQC
    from cpg_workflows.workflow import get_workflow
    from cpg_workflows.stages.seqr_loader import MtToEs

    get_workflow().run(stages=[MtToEs, GvcfMultiQC, CramMultiQC, JointVcfQC])

    assert (
        get_batch().job_by_tool['gatk HaplotypeCaller']['job_n']
        == len(get_cohort().get_samples()) * 50
    )
    assert get_batch().job_by_tool['picard MergeVcfs']['job_n'] == len(
        get_cohort().get_samples()
    )
    assert get_batch().job_by_tool['gatk ReblockGVCF']['job_n'] == len(
        get_cohort().get_samples()
    )
    assert (
        get_batch().job_by_tool['picard CollectVariantCallingMetrics']['job_n']
        == len(get_cohort().get_samples()) + 1
    )
    assert get_batch().job_by_tool['gatk GenomicsDBImport']['job_n'] == 50
