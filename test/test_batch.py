"""
Test initializing Batch object.
"""

import hail as hl
import toml

from cpg_utils import to_path, Path
from cpg_utils.config import set_config_paths, update_dict
from cpg_utils.hail_batch import dataset_path, command
from cpg_workflows.utils import timestamp
from cpg_workflows.batch import get_batch

tmp_dir_path = to_path(__file__).parent / 'results' / timestamp()
tmp_dir_path = tmp_dir_path.absolute()
tmp_dir_path.mkdir(parents=True, exist_ok=True)

DEFAULT_CONF = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
sequencing_type = 'genome'

check_inputs = false
check_intermediates = false
check_expected_outputs = false

[storage.default]
default = '{str(tmp_dir_path)}'

[storage.fewgenomes]
default = '{str(tmp_dir_path)}'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
"""


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    if extra_conf:
        update_dict(d, extra_conf)
    config_path = dir_path / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(config_path)])


def test_batch_job():
    """
    Test creating a job and running a batch.
    """
    _set_config(tmp_dir_path)
    b = get_batch('Test batch job')
    j1 = b.new_job('Jo b1')
    text = 'success'
    cmd = f"""\
    echo {text} > {j1.output}
    """
    j1.command(command(cmd))
    output1_path = dataset_path('output1.txt')
    b.write_output(j1.output, str(output1_path))

    j2 = b.new_job('Job 2')
    j2.command(f'touch {j2.output}')
    j2.command(f'cat {b.read_input(output1_path)} >> {j2.output}')
    j2.depends_on(j1)
    output2_path = dataset_path('output2.txt')
    b.write_output(j2.output, str(output2_path))

    b.run()
    with to_path(output2_path).open() as fh:
        assert fh.read().strip() == text
