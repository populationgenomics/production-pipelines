"""
Providers that can be overriden in specific infrastructutres.
"""

from enum import Enum


class InputProviderType(Enum):
    """Available types of input prpviders"""

    SMDB = 'smdb'
    CSV = 'csv'


class StoragePolicyType(Enum):
    """Available storage policies"""

    CPG = 'cpg'


class StatusReporterType(Enum):
    """Available types of status reporters"""

    SMDB = 'smdb'
