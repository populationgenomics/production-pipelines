"""
Types specific to Cloud storage. 
"""
from enum import Enum
from abc import ABC, abstractmethod

from cloudpathlib import CloudPath


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
    ) -> CloudPath:
        """
        Bucket to write results:
        @param dataset: dataset/stack name
        @param namespace: namespace (test or main)
        @param suffix: (optional) suffix to append to the bucket name
        @param version: (optional) pipeline version
        @param sample: (optional) sample name
        """

    def get_analysis_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> CloudPath:
        """
        Bucket for analysis results.
        """
        return self.get_bucket(
            dataset=dataset, 
            namespace=namespace, 
            version=version, 
            sample=sample, 
            suffix='analysis'
        )

    def get_tmp_bucket(
        self, 
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> CloudPath:
        """
        Bucket for temporary files.
        """
        return self.get_bucket(
            dataset=dataset, 
            namespace=namespace, 
            version=version, 
            sample=sample, 
            suffix='tmp'
        )

    def get_web_bucket(
        self, 
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> CloudPath:
        """
        Bucket shared with an HTTP server.
        """
        return self.get_bucket(
            dataset=dataset, 
            namespace=namespace, 
            version=version, 
            sample=sample, 
            suffix='web'
        )

    # noinspection PyMethodMayBeStatic
    def get_web_url(
        self, 
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> str | None:
        """
        URL corrsponding to the WEB bucket.
        """
        return None
