"""
CPG implementation of StorageProvider:
https://github.com/populationgenomics/team-docs/tree/main/storage_policies
across Google Cloud Storage and Azure Blob Storage.
"""
from cloudpathlib import CloudPath

from cpg_pipes.storage import Path

from cpg_pipes.storage import StorageProvider, Cloud, Namespace


class CPGStorageProvider(StorageProvider):
    """
    CPG storage policy implementation of the StorageProvider
    """
    def __init__(self, cloud: Cloud = Cloud.GS):
        super().__init__(cloud)
        self.prefix = 'cpg'

    def get_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Bucket name is constructed according to the storage policy:
        https://github.com/populationgenomics/team-docs/tree/main/storage_policies
        """
        path = CloudPath(
            f'{self.cloud.value}://'
            f'{self.prefix}-{dataset}-{namespace.value}'
        )
        if suffix:
            path = CloudPath(f'{path}-{suffix}')
        if version:
            path = path / version
        if sample:
            path = path / sample
        return path

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
        url = (
            f'https://{namespace.value}-web.populationgenomics.org.au/'
            f'{dataset}'
        )
        if version:
            url += f'/{version}'
        if sample:
            url += f'/{sample}'
        return url
