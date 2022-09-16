"""
Test Hail Query functions.
"""
import shutil
import pytest
import toml
from cpg_utils import to_path, Path
from cpg_utils.config import update_dict
from cpg_utils.workflows.batch import get_batch
from cpg_utils.workflows.filetypes import BamPath, FastqPair, FastqPairs
from cpg_utils.workflows.targets import Cohort
from cpg_utils.workflows.utils import timestamp
from cpg_utils.workflows.workflow import get_workflow
from cpg_utils.config import set_config_paths
from stages.seqr_loader import MtToEs

DEFAULT_CONF = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
check_inputs = false
check_intermediates = false
check_expected_outputs = false
assume_outputs_exist = true
sequencing_type = 'genome'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = true
"""


@pytest.fixture()
def tmp_dir() -> Path:
    dir_path = to_path('results') / timestamp()
    dir_path.mkdir(parents=True, exist_ok=True)
    yield dir_path
    shutil.rmtree(dir_path)


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    if extra_conf:
        update_dict(d, extra_conf)
    config_path = dir_path / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(config_path)])


def test_workflow_dry(mocker, tmp_dir: Path):
    """
    Run entire seqr-loader in a dry mode.
    """
    _set_config(
        tmp_dir,
        extra_conf={
            'workflow': {
                'dataset': 'test-analysis-dataset',
                'skip_stages': ['Align'],
                'last_stages': ['AnnotateDataset'],
            },
            'hail': {
                'dry_run': True,
                # 'backend': 'local',  # j.cloudfuse doesn't work with local backend
            },
        },
    )

    cohort = Cohort()
    ds = cohort.create_dataset('test-input-dataset')
    ds.add_sample('CPG01', 'SAMPLE1', alignment_input_by_seq_type={
        'genome': BamPath('gs://test-input-dataset-upload/sample1.bam')
    })
    ds.add_sample('CPG02', 'SAMPLE2', alignment_input_by_seq_type={
        'genome': FastqPairs([
            FastqPair(
                'gs://test-input-dataset-upload/sample2_L1_R1.fq.gz', 
                'gs://test-input-dataset-upload/sample2_L1_R2.fq.gz'
            ),
            FastqPair(
                'gs://test-input-dataset-upload/sample2_L2_R1.fq.gz', 
                'gs://test-input-dataset-upload/sample2_L2_R2.fq.gz'
            ),
        ])
    })
    
    def mock_exists(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> bool:
        return False

    # mocker.patch('cpg_utils.workflows.utils.exists', mock_exists)
    mocker.patch('cloudpathlib.cloudpath.CloudPath.exists', mock_exists)

    def mock_get_cohort(  # pylint: disable=unused-argument
        *args, **kwargs
    ) -> Cohort:
        return cohort

    mocker.patch('cpg_utils.workflows.inputs.create_cohort', mock_get_cohort)

    get_workflow().run(stages=[MtToEs])
    
    assert (
        get_batch().job_by_tool['gatk HaplotypeCaller']['job_n']
        == len(cohort.get_samples()) * 50
    )
    assert get_batch().job_by_tool['picard MergeVcfs']['job_n'] == len(cohort.get_samples())
    assert get_batch().job_by_tool['gatk ReblockGVCF']['job_n'] == len(cohort.get_samples())
    assert (
        get_batch().job_by_tool['picard CollectVariantCallingMetrics']['job_n']
        == len(cohort.get_samples()) + 1
    )
    assert get_batch().job_by_tool['gatk GenomicsDBImport']['job_n'] == 2
