"""
Types specific to Cloud storage. 

The default assumption is the CPG storage policy:
https://github.com/populationgenomics/team-docs/tree/main/storage_policies
across Google Cloud Storage and Azure Blob Storage.
"""
from enum import Enum


class Namespace(Enum):
    """
    CPG storage namespace. See for more details on storage policies:
    https://github.com/populationgenomics/team-docs/tree/main/storage_policies#main-vs-test
    """
    MAIN = 'main'  # production runs: read from main, write to main
    TEST = 'test'  # runs that make test data: read from test, write to test
    TMP = 'tmp'    # read from test, write to tmp


class StorageProvider(Enum):
    """
    Cloud infrastructure provider.
    """
    GS = 'gs'
    AZ = 'az'
