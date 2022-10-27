"""
Test Hail Query functions.
"""
import shutil
import pytest
import toml
from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths, update_dict
from cpg_workflows.utils import timestamp
from cpg_workflows.python_scripts import check_multiqc

DEFAULT_CONF = """
[workflow]
sequencing_type = 'genome'

[qc_thresholds.genome.min]
"MEDIAN_COVERAGE" = 10
"PCT_PF_READS_ALIGNED" = 80
[qc_thresholds.genome.max]
"FREEMIX" = 0.04
"PERCENT_DUPLICATION" = 25
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


def test_check_multiqc(caplog, tmp_dir: Path):
    _set_config(tmp_dir)
    data_dir = to_path(__file__).parent / 'data' / 'check_multiqc'
    check_multiqc.run(str(data_dir / 'validation_multiqc.json'))
    for expected_line in ['⭕ CPG243717|NA12878_KCCG: MEDIAN_COVERAGE=8.00<10.00']:
        matching_lines = [expected_line in msg for msg in caplog.messages]
        assert any(matching_lines)
