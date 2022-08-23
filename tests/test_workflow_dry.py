"""
Test Hail Query functions.
"""
import shutil

import hail as hl
import pytest
import toml
from cpg_utils import to_path, Path
from cpg_utils.config import get_config, set_config_paths, update_dict
from cpg_utils.workflows.batch import get_batch
from cpg_utils.workflows.utils import timestamp
from cpg_utils.hail_batch import dataset_path, init_batch


DEFAULT_CONF = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
check_inputs = false
check_intermediates = false
check_expected_outputs = false
jc_intervals_num = 2
vqsr_intervals_num = 2
vep_intervals_num = 2
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


def _local_backend_config(dir_path: Path):
    _set_config(dir_path, {
        'workflow': {
            'path_scheme': 'local', 
            'local_dir': str(dir_path),
        },
        'hail': {'backend': 'local'},
    })


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    if extra_conf:
        update_dict(d, extra_conf)
    config_path = dir_path / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(config_path)])


def test_workflow_dry(tmp_dir: Path):
    """
    Run entire seqr-loader in a dry mode.
    """
    intervals_num = 2
    sample_ids = ['CPG56564', 'CPG56572']
    _set_config(tmp_dir, extra_conf={
        'workflow': {
            'stages': 'MtToEs',
            'only_samples': sample_ids,
            'skip_stages': ['Align'],
            'hc_intervals_num': intervals_num,
            'jc_intervals_num': intervals_num,
        },
        'hail': {
            'dry_run': True,
        }
    })
    import main
    main.main()
    assert get_batch().job_by_tool['gatk_HaplotypeCaller']['job_n'] == len(sample_ids) * intervals_num
    assert get_batch().job_by_tool['picard_MergeVcfs']['job_n'] == len(sample_ids)
    assert get_batch().job_by_tool['gatk_ReblockGVCF']['job_n'] == len(sample_ids)
    assert get_batch().job_by_tool['picard_CollectVariantCallingMetrics']['job_n'] == len(sample_ids) + 1
    assert get_batch().job_by_tool['gatk_GenomicsDBImport']['job_n'] == intervals_num
