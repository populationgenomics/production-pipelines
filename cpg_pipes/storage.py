"""
Types specific to Cloud storage. 

The default assumption is the CPG storage policy:
https://github.com/populationgenomics/team-docs/tree/main/storage_policies
across Google Cloud Storage and Azure Blob Storage.
"""
from enum import Enum
from abc import ABC, abstractmethod

from cloudpathlib import CloudPath


class Namespace(Enum):
    """
    CPG storage namespace. See for more details on storage policies:
    https://github.com/populationgenomics/team-docs/tree/main/storage_policies#main-vs-test
    """
    MAIN = 'main'
    TEST = 'test'


class Cloud(Enum):
    """
    Cloud storage protocol.
    """
    GS = 'gs'
    AZ = 'az'


class StorageProvider(ABC):
    """
    Abstract class for implementing storage bucket policy.
    Onty get_bucket() method is required, however other methods
    are available to override as well.
    """
    def __init__(self, protocol: Cloud):
        self.protocol = protocol

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
    
    
class CPGStorageProvider(StorageProvider):
    """
    CPG storage policy implementation of the StorageProvider
    """
    def __init__(self, prefix: Cloud):
        super().__init__(prefix)
        self.prefix = 'cpg'
        
    def get_bucket(
        self, 
        dataset: str,
        namespace: Namespace,
        suffix: str = None,
        version: str | None = None,
        sample_name: str = None,
    ) -> CloudPath:
        """
        Bucket name is constructed according to the storage policy:
        https://github.com/populationgenomics/team-docs/tree/main/storage_policies
        """
        path = CloudPath(
            f'{self.protocol.value}://'
            f'{self.prefix}-{dataset}-{namespace.value}'
        )
        if suffix:
            path = CloudPath(f'{path}-{suffix}')
        if version:
            path = path / version
        if sample_name:
            path = path / sample_name
        return path

    # noinspection PyMethodMayBeStatic
    def get_web_url(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample_name: str = None,
    ) -> str | None:
        """
        URL corrsponding to the WEB bucket.
        """
        url = (
            f'https://{namespace.value}-web.populationgenomics.org.au/'
            f'{dataset}'
        )
        if version:
            url += f'/{version}'
        if sample_name:
            url += f'/{sample_name}'
        return url
