"""
CPG implementation of StorageProvider, according to the storage policies:
https://github.com/populationgenomics/team-docs/tree/main/storage_policies
Works only on Google Cloud Storage for now.
"""

from cloudpathlib import CloudPath

from ..storage import StorageProvider, Cloud, Namespace, Path, to_path


class CpgStorageProvider(StorageProvider):
    """
    CPG storage policy implementation of the StorageProvider
    """

    def __init__(self, cloud: Cloud = Cloud.GS):
        self.gc_prefix = 'cpg'
        super().__init__(cloud)

    def _dataset_base(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str = None,
    ) -> Path:
        """
        Base path for a dataset.
        """
        container = f'{dataset}-{namespace.value}'
        if suffix:
            container = f'{container}-{suffix}'
        return CloudPath(f'{self.cloud.value}://{self.gc_prefix}-{container}')

    def get_base(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Bucket name is constructed according to the CPG storage policy:
        https://github.com/populationgenomics/team-docs/tree/main/storage_policies
        """
        path = self._dataset_base(dataset, namespace, suffix)
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
        URL corresponding to the WEB bucket.
        """
        url = f'https://{namespace.value}-web.populationgenomics.org.au/{dataset}'
        if version:
            url += f'/{version}'
        if sample:
            url += f'/{sample}'
        return url

    def get_ref_base(self) -> Path:
        """
        Prefix for reference data.
        """
        return to_path(f'{self.cloud.value}://cpg-reference')
