"""
Test entire seqr-loader in a dry mode.
"""
import os

import toml
from cpg_utils import to_path, Path
from cpg_utils.config import update_dict
from cpg_utils.workflows.batch import get_batch
from cpg_utils.workflows.filetypes import BamPath, FastqPair, FastqPairs
from cpg_utils.workflows.inputs import get_cohort
from cpg_utils.workflows.targets import Cohort
from cpg_utils.workflows.utils import timestamp
from cpg_utils.workflows.workflow import get_workflow
from cpg_utils.config import set_config_paths
from cpg_workflows.stages.seqr_loader import MtToEs
from pytest_mock import MockFixture


def _set_config(results_prefix: Path, extra_conf: dict | None = None):
    with (
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'seqr_loader.toml'
    ).open() as f:
        d = toml.load(f)
    update_dict(
        d,
        {
            'workflow': {
                'dataset_gcp_project': 'test-analysis-dataset-1234',
                'dataset': 'test-analysis-dataset',
                'access_level': 'test',
                'sequencing_type': 'genome',
                'driver_image': '<stub>',
                'skip_stages': ['Align'],
                'check_inputs': False,
                'check_intermediates': False,
                'check_expected_outputs': False,
                'path_scheme': 'local',
                'status_reporter': None,
            },
            'hail': {
                'billing_project': 'test-analysis-dataset',
                'delete_scratch_on_exit': True,
                'dry_run': True,
                'backend': 'local',
            },
        },
    )
    if extra_conf:
        update_dict(d, extra_conf)
    with (out_path := results_prefix / 'config.toml').open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(out_path)])


def test_seqr_loader_dry(mocker: MockFixture):
    """
    Test entire seqr-loader in a dry mode.
    """
    results_prefix = (
        to_path(__file__).parent / 'results' / os.getenv('TEST_TIMESTAMP', timestamp())
    ).absolute()
    results_prefix.mkdir(parents=True, exist_ok=True)

    _set_config(
        results_prefix=results_prefix,
        extra_conf={
            'local_dir': str(results_prefix),
        },
    )

    def mock_create_cohort(*args, **kwargs) -> Cohort:
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

    def mock_exists(*args, **kwargs) -> bool:
        return False

    def do_nothing(*args, **kwargs):
        return None

    mocker.patch('cpg_utils.workflows.inputs.create_cohort', mock_create_cohort)
    # functions like get_intervals checks file existence
    mocker.patch('cloudpathlib.cloudpath.CloudPath.exists', mock_exists)
    # cloudfuse (used in Vep) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.cloudfuse', do_nothing)
    # always_run (used in MtToEs -> hail_dataproc_job) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.always_run', do_nothing)
    # can't access secrets from CI environment
    mocker.patch('stages.seqr_loader.es_password', lambda: 'test-password')

    get_workflow().run(stages=[MtToEs])

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
