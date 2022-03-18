import coloredlogs
from os.path import dirname, abspath

coloredlogs.install(
    level='DEBUG', 
    fmt='%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
)


def get_package_path():
    """
    @return: local install path of the package
    """
    return dirname(abspath(__file__))
