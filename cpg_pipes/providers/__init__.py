"""
Providers that can be overridden in specific infrastructures.
"""

from enum import Enum


class InputProviderType(Enum):
    """Available types of input providers"""

    CPG = 'cpg'
    CSV = 'csv'


class StatusReporterType(Enum):
    """Available types of status reporters"""

    CPG = 'cpg'
