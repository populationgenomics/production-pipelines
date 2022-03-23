"""
Providers that can be overriden in specific infrastructutres.
"""

from enum import Enum

from .storage import Namespace, Cloud


class InputProviderType(Enum):
    """Available types of input prpviders"""
    SMDB = 'smdb'
    CSV = 'csv'
    NONE = 'none'


class StoragePolicy(Enum):
    """Available storage policies"""
    CPG = 'cpg'


class StatusReporterType(Enum):
    """Available types of status reporters"""
    SMDB = 'smdb'
    NONE = 'none'
