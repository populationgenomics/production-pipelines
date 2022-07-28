import coloredlogs
import importlib.metadata

from cpg_utils.hail_batch import Path, to_path, Namespace

__all__ = [
    'Path', 
    'to_path', 
    'Namespace', 
    'get_package_name', 
    'get_package_path', 
    'get_version'
]


coloredlogs.install(
    level='DEBUG', fmt='%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
)


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
