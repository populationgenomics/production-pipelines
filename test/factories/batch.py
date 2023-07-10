from pathlib import Path

from hailtop.batch import LocalBackend

from cpg_workflows.batch import Batch


def create_local_batch(tmp_path: Path | str) -> Batch:
    """
    FIXME: __init__ looks in config for hail `dry_run` so this function needs to run
    after config is set up via `set_config_paths` to avoid a `KeyError`.
    """
    return Batch(name='test-align_job', backend=LocalBackend(tmp_dir=str(tmp_path)))
