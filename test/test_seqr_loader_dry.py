"""
Test seqr-loader workflow.
"""

from pathlib import Path
from unittest.mock import mock_open

from pytest_mock import MockFixture

from cpg_utils import to_path

from . import set_config

TOML = """
[workflow]
dataset_gcp_project = "test-analysis-dataset-1234"
dataset = "test-analysis-dataset"
access_level = "test"
sequencing_type = "genome"
driver_image = "<stub>"
skip_stages = ["Align"]
check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = "local"
local_dir = "{directory}"

[hail]
billing_project = "test-analysis-dataset"
delete_scratch_on_exit = true
dry_run = true
backend = "local"

[images]
cpg_workflows = "stub"
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
vep_110 = "stub"
vep_105 = "stub"
verifybamid = "stub"

[references]
genome_build = "GRCh38"
vep_mount = "stub"
liftover_38_to_37 = "stub"
somalier_sites = "stub"
seqr_combined_reference_data = "stub"
seqr_clinvar = "stub"
[references.hg38_telomeres_and_centromeres_intervals]
interval_list = "stub"

[storage.default]
default = "{directory}"
web = "{directory}-web"
analysis = "{directory}-analysis"
tmp = "://{directory}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.test-input-dataset]
default = "{directory}"
web = "{directory}-web"
analysis = "{directory}-analysis"
tmp = "{directory}-test-tmp"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.test-analysis-dataset]
default = "{directory}"
web = "{directory}-web"
analysis = "{directory}-analysis"
tmp = "{directory}-test-tmp"
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

DEFAULT_CONFIG = Path(to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml')
SEQR_LOADER_CONFIG = Path(to_path(__file__).parent.parent / 'configs' / 'defaults' / 'seqr_loader.toml')


def _mock_cohort():
    from cpg_workflows.filetypes import BamPath, FastqPair, FastqPairs
    from cpg_workflows.targets import Cohort

    cohort = Cohort()
    ds = cohort.create_dataset('test-input-dataset')
    ds.add_sequencing_group(
        'CPGAA',
        'SAMPLE1',
        alignment_input_by_seq_type={'genome': BamPath('gs://test-input-dataset-upload/sample1.bam')},
    )
    ds.add_sequencing_group(
        'CPGBB',
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
                ],
            ),
        },
    )
    return cohort


def selective_mock_open(*args, **kwargs):
    if str(args[0]).endswith('.toml'):
        # Don't mock calls to load a config file
        return open(*args, **kwargs)
    else:
        return mock_open(read_data='<stub>')(*args, **kwargs)


def test_seqr_loader_dry(mocker: MockFixture, tmp_path):
    """
    Test entire seqr-loader in a dry mode.
    """
    conf = TOML.format(directory=str(tmp_path))
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG, SEQR_LOADER_CONFIG],
    )

    mocker.patch('cpg_workflows.inputs.deprecated_create_cohort', _mock_cohort)

    def do_nothing(*args, **kwargs):
        return None

    def mock_create_analysis(*args, **kwargs) -> int:
        return 1

    mocker.patch('pathlib.Path.open', selective_mock_open)
    # cloudfuse (used in Vep) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.cloudfuse', do_nothing)
    # always_run (used in MtToEs -> hail_dataproc_job) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.always_run', do_nothing)
    # can't access secrets from CI environment
    mocker.patch('cpg_workflows.stages.seqr_loader.es_password', lambda: 'test-password')
    mocker.patch(
        'metamist.apis.AnalysisApi.create_analysis',
        mock_create_analysis,
    )
    mocker.patch('metamist.apis.AnalysisApi.update_analysis', do_nothing)

    from cpg_utils.hail_batch import get_batch
    from cpg_workflows.inputs import get_multicohort
    from cpg_workflows.stages.cram_qc import CramMultiQC
    from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
    from cpg_workflows.stages.joint_genotyping_qc import JointVcfQC
    from cpg_workflows.stages.seqr_loader import MtToEs
    from cpg_workflows.workflow import get_workflow

    get_workflow().run(stages=[MtToEs, GvcfMultiQC, CramMultiQC, JointVcfQC])

    assert (
        get_batch().job_by_tool['gatk HaplotypeCaller']['job_n'] == len(get_multicohort().get_sequencing_groups()) * 50
    )
    assert get_batch().job_by_tool['picard MergeVcfs']['job_n'] == len(get_multicohort().get_sequencing_groups())
    assert get_batch().job_by_tool['gatk ReblockGVCF']['job_n'] == len(get_multicohort().get_sequencing_groups())
    assert (
        get_batch().job_by_tool['picard CollectVariantCallingMetrics']['job_n']
        == len(get_multicohort().get_sequencing_groups()) + 1
    )
    assert get_batch().job_by_tool['gatk GenomicsDBImport']['job_n'] == 50
