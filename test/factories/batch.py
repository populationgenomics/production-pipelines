from pathlib import Path

from hailtop.batch import LocalBackend

from cpg_workflows.batch import Batch


def create_local_batch(tmp_path: Path | str) -> Batch:
    return Batch(name="test-align_job", backend=LocalBackend(tmp_dir=str(tmp_path)))
