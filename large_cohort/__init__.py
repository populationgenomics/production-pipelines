import importlib.metadata
import logging

from cpg_utils import Path, to_path

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


def get_package_name() -> str:
    """
    Get name of the package.
    """
    return __name__.split('.', 1)[0]


def get_package_path() -> Path:
    """
    Get local install path of the package.
    """
    return to_path(__file__).parent.absolute()


def get_version() -> str:
    """
    Get package version.
    """
    return importlib.metadata.version(get_package_name())
