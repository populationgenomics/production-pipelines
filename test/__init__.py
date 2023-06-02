import logging
from pathlib import Path
from typing import Any

import toml
from cpg_utils.config import set_config_paths

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


def update_dict(d1: dict, d2: dict) -> None:
    """
    Helper functon to merge one dict into another. We don't import `update_dict` from
    `cpg_utils`, as it would break mocking in some tests.
    """
    for k, v2 in d2.items():
        v1 = d1.get(k)
        if isinstance(v1, dict) and isinstance(v2, dict):
            update_dict(v1, v2)
        else:
            d1[k] = v2


def set_config(config: str | dict[str, Any], path: Path, merge_with: list[Path] = []):
    with path.open('w') as f:
        if isinstance(config, dict):
            toml.dump(config, f)
        else:
            f.write(config)

        f.flush()

    set_config_paths([*[str(s) for s in merge_with], str(path)])
