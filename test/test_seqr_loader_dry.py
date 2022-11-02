"""
Test entire seqr-loader in a dry mode.
"""
import os
import string
import time
from random import choices
from unittest.mock import mock_open

import toml
from cpg_utils import to_path, Path
from pytest_mock import MockFixture


def update_dict(d1: dict, d2: dict) -> None:
    """Updates the d1 dict with the values from the d2 dict recursively in-place."""
    for k, v2 in d2.items():
        v1 = d1.get(k)
        if isinstance(v1, dict) and isinstance(v2, dict):
            update_dict(v1, v2)
        else:
            d1[k] = v2


def timestamp(rand_suffix_len: int = 5) -> str:
    """
    Generate a timestamp string. If `rand_suffix_len` is set, adds a short random
    string of this length for uniqueness.
    """
    result = time.strftime('%Y_%m%d_%H%M')
    if rand_suffix_len:
        rand_bit = ''.join(
            choices(string.ascii_uppercase + string.digits, k=rand_suffix_len)
        )
        result += f'_{rand_bit}'
    return result


def _make_config(results_prefix: Path) -> dict:
    d: dict = {}
    for fp in (
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'workflows.toml',
        to_path(__file__).parent.parent / 'configs' / 'defaults' / 'seqr_loader.toml',
    ):
        with fp.open():
            update_dict(d, toml.load(fp))

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
                'local_dir': str(results_prefix),
            },
            'hail': {
                'billing_project': 'test-analysis-dataset',
                'delete_scratch_on_exit': True,
                'dry_run': True,
                'backend': 'local',
            },
        },
    )
    return d


def test_seqr_loader_dry(mocker: MockFixture):
    """
    Test entire seqr-loader in a dry mode.
    """
    results_prefix = (
        to_path(__file__).parent / 'results' / os.getenv('TEST_TIMESTAMP', timestamp())
    ).absolute()
    results_prefix.mkdir(parents=True, exist_ok=True)
    conf = _make_config(results_prefix)

    def mock_create_cohort(*args, **kwargs):
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

    def mock_exists(*args, **kwargs) -> bool:
        return False

    def do_nothing(*args, **kwargs):
        return None

    def mock_create_new_analysis(*args, **kwargs) -> int:
        return 1

    mocker.patch('pathlib.Path.open', mock_open(read_data='<stub>'))
    mocker.patch('cpg_utils.config.get_config', lambda: conf)
    mocker.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)
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
    from cpg_workflows.workflow import get_workflow
    from cpg_workflows.stages.seqr_loader import MtToEs

    get_workflow().run(stages=[MtToEs, GvcfMultiQC, CramMultiQC])

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
