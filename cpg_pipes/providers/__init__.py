"""
Providers that can be overriden in specific infrastructutres.
"""

from enum import Enum


class InputProviderType(Enum):
    """Available types of input prpviders"""

    SMDB = 'smdb'
    CSV = 'csv'


class StatusReporterType(Enum):
    """Available types of status reporters"""

    SMDB = 'smdb'
