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


class Cloud(Enum):
    """
    Cloud storage provider and correponding protocol prefix.
    """

    GS = 'gs'
    AZ = 'az'


class StorageProvider(ABC):
    """
    Abstract class for implementing storage bucket policy.
    Onty get_bucket() method is required, however other methods
    are available to override as well.
    """

    def __init__(self, cloud: Cloud):
        self.cloud = cloud

    @abstractmethod
    def get_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str | None = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Bucket to write primary results.
        @param dataset: dataset/stack name
        @param namespace: namespace (test or main)
        @param suffix: (optional) suffix to append to the bucket name
        @param version: (optional) pipeline version
        @param sample: (optional) sample name
        """

    def get_tmp_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Bucket to write temporary files.
        """
        return self.get_bucket(
            dataset=dataset,
            namespace=namespace,
            version=version,
            sample=sample,
            suffix='tmp',
        )

    def get_web_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Bucket shared with an HTTP server.
        """
        return self.get_bucket(
            dataset=dataset,
            namespace=namespace,
            version=version,
            sample=sample,
            suffix='web',
        )

    @abstractmethod
    def get_web_url(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> str | None:
        """
        URL corresponding to the WEB bucket.
        """

    @abstractmethod
    def get_ref_bucket(self) -> Path:
        """
        Prefix for reference data.
        """
