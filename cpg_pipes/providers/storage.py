"""
Abstract storage provider.
"""
import pathlib
from enum import Enum
from abc import ABC, abstractmethod
from typing import Union

from cloudpathlib import CloudPath
from cloudpathlib.anypath import to_anypath


# Path can be either a cloud URL or a local posix file path.
Path = Union[CloudPath, pathlib.Path]

# Using convenience method from cloudpathlib to parse a path string.
to_path = to_anypath


class Namespace(Enum):
    """
    Storage namespace.
    https://github.com/populationgenomics/team-docs/tree/main/storage_policies#main-vs-test
    """

    MAIN = 'main'
    TEST = 'test'


class StorageProvider(ABC):
    """
    Abstract class for implementing storage path policy.
    Onty path() method is required, however other methods are available to 
    override as well.
    """

    def __init__(self):
        pass

    @abstractmethod
    def path(
        self,
        dataset: str,
        namespace: Namespace,
        category: str | None = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Path prefix for primary results.
        @param dataset: dataset name
        @param namespace: namespace (test or main)
        @param category: (optional) category to append to the bucket name
        @param version: (optional) pipeline version
        @param sample: (optional) sample name
        """

    @abstractmethod
    def web_url(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> str | None:
        """
        URL corresponding to path(category='web', ...), assuming other parameters
        are the same.
        """
