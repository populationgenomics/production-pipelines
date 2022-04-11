"""
Providers that can be overridden in specific infrastructures.
"""

from enum import Enum


class InputProviderType(Enum):
    """Available types of input providers"""

    SMDB = 'smdb'
    CSV = 'csv'


class StoragePolicyType(Enum):
    """Available storage policies"""

    CPG = 'cpg'


class StatusReporterType(Enum):
    """Available types of status reporters"""

    SMDB = 'smdb'
