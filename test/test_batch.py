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
path_scheme = 'local'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
"""


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    d['workflow']['local_dir'] = str(dir_path)
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


def test_batch_python_job():
    """
    Testing calling a python job.
    """
    _set_config(tmp_dir_path)

    b = get_batch('Test batch python job')
    j = b.new_python_job('Test python job')

    input_tsv_path = to_path(dataset_path('input.tsv'))
    input_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with input_tsv_path.open('w') as f:
        f.write('col1\tcol2\n1\t2')

    def query_fn(tsv_path: str, out_ht_path: str):
        ht = hl.import_table(tsv_path, types={'col1': hl.tint, 'col2': hl.tint})
        ht.show()
        ht = ht.annotate(col3=ht.col1 + ht.col2)
        ht.write(out_ht_path, overwrite=True)

    output_ht_path = dataset_path('output.ht')
    j.call(query_fn, str(input_tsv_path), output_ht_path)
    b.run()

    hl.init_local(log=dataset_path('hail-log.txt'))
    result = hl.read_table(str(output_ht_path)).col3.collect()[0]
    assert result == 3, result
