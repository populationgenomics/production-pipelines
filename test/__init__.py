import logging
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import toml

from cpg_utils import Path as AnyPath
from cpg_utils.config import set_config_paths

logging.basicConfig()
logging.getLogger().setLevel(logging.WARN)

# silence other loggers
logging.getLogger('pyspark').setLevel(logging.ERROR)
logging.getLogger('py4j').setLevel(logging.ERROR)


@runtime_checkable
class IDictRepresentable(Protocol):
    def as_dict(self) -> dict[str, Any]: ...


class TomlAnyPathEncoder(toml.TomlEncoder):
    """
    Support for CPG path objects in TOML, which might be regular a pathlib path or
    a cloud path.
    """

    def dump_value(self, v):
        if isinstance(v, AnyPath):
            v = str(v)
        return super().dump_value(v)


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


def set_config(
    config: str | dict[str, Any] | IDictRepresentable,
    path: Path,
    merge_with: list[Path] | None = None,
) -> None:
    """
    Writes your config to `path` and sets the `CPG_CONFIG_PATH` environment variable
    so that the `cpg_utils` module will be able to read your config. If `merge_with` is
    provided, the config will be merged with the configs at the given paths. Merging
    happens right to left, so that values in the right config will override values
    in the left config.

    Args:
        config (str | dict[str, Any] | IDictRepresentable):
            A valid TOML string, a dictionary to be converted to TOML, or an object
            which implements the `IDictRepresentable` protocol.

        path (Path):
            Path to write the config to.

        merge_with (list[Path] | None, optional):
            A list of paths to merge with the config. Merging happens right to left,
            so that values in the right config will override values in the left config.
            Defaults to `None`.
    """
    with path.open('w') as f:
        if isinstance(config, dict):
            toml.dump(config, f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, IDictRepresentable):
            toml.dump(config.as_dict(), f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, str):
            f.write(config)
        else:
            raise TypeError(f'Expected config to be a string, dict, or IDictRepresentable, butgot {type(config)}')

        f.flush()

    return set_config_paths([*[str(s) for s in (merge_with or [])], str(path)])
