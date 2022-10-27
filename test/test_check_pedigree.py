"""
Test Hail Query functions.
"""

import shutil
import pytest
import toml
from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths, update_dict
from cpg_workflows.utils import timestamp
from cpg_workflows.python_scripts import check_pedigree

DEFAULT_CONF = """
[workflow]
sequencing_type = 'genome'
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


def test_check_pedigree(caplog, tmp_dir: Path):
    _set_config(tmp_dir)
    data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'
    check_pedigree.run(
        somalier_samples_fpath=str(data_dir / 'somalier-samples.tsv'),
        somalier_pairs_fpath=str(data_dir / 'somalier-pairs.tsv'),
        expected_ped_fpath=str(data_dir / 'samples.ped'),
    )
    for expected_line in [
        '4/51 PED samples with mismatching sex',
        'Sex inferred for 51/51 samples, matching for 47 samples.',
        'Found 4 sample pair(s) that are provided as unrelated, are inferred as related:',
        'Found 7 sample pair(s) that are provided as related, but inferred as unrelated:',
    ]:
        matching_lines = [expected_line in msg for msg in caplog.messages]
        assert any(matching_lines)
