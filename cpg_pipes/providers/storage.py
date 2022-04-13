"""
Abstract storage provider.
"""
import os
import pathlib
from enum import Enum
from abc import ABC, abstractmethod
from typing import Union

from cloudpathlib import CloudPath
from cloudpathlib.anypath import to_anypath
import hail_az  # importing to register hail-az prefix.


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
    Cloud storage provider and corresponding protocol prefix.
    """

    GS = 'gs'
    HAIL_AZ = 'hail-az'


class StorageProvider(ABC):
    """
    Abstract class for implementing storage path policy.
    Onty get_base() method is required, however other methods are available to 
    override as well.
    """

    def __init__(
        self, 
        cloud: Cloud, 
        az_account: str | None = None,
        az_container: str | None = None,
    ):
        self.cloud = cloud
        self.az_account = az_account or os.environ.get('STORAGE_ACCOUNT_NAME')
        self.az_container = az_container or os.environ.get('STORAGE_CONTAINER_NAME')

    @abstractmethod
    def get_base(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str | None = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Base path for primary results.
        @param dataset: dataset/stack name
        @param namespace: namespace (test or main)
        @param suffix: (optional) suffix to append to the bucket name
        @param version: (optional) pipeline version
        @param sample: (optional) sample name
        """

    def get_tmp_base(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Base path for temporary files.
        """
        return self.get_base(
            dataset=dataset,
            namespace=namespace,
            version=version,
            sample=sample,
            suffix='tmp',
        )

    def get_web_base(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Path base corresponding to the HTTP server. Used to expose files in web:
        `get_web_url()` with same parameters should return a HTTP URL matching this 
        object. Inspired by the CPG storage policy: 
        https://github.com/populationgenomics/team-docs/blob/main/storage_policies/README.md#web-gscpg-dataset-maintest-web
        """
        return self.get_base(
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
    def get_ref_base(self) -> Path:
        """
        Prefix for reference data.
        """
