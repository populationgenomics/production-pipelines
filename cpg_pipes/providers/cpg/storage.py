"""
CPG storage policy implementation of the StorageProvider
Using cpg_utils to construct paths, according to the CPG storage policy:
https://github.com/populationgenomics/team-docs/tree/main/storage_policies
"""

from cpg_utils import hail_batch

from ..storage import StorageProvider, Namespace, Path, to_path


class CpgStorageProvider(StorageProvider):
    """
    CPG storage policy implementation of the StorageProvider
    Using cpg_utils to construct paths, according to the CPG storage policy:
    https://github.com/populationgenomics/team-docs/tree/main/storage_policies
    """

    # noinspection PyMethodMayBeStatic
    def path(
        self,
        dataset: str,
        namespace: Namespace,
        category: str = None,
        version: str | None = None,
        sample: str = None,
    ) -> Path:
        """
        Path prefix for primary results.
        """
        path = to_path(hail_batch.dataset_path(
            suffix='', 
            category=category,
            dataset=dataset,
            access_level=namespace.value,
        ))
        if version:
            path = path / version
        if sample:
            path = path / sample
        return path

    # noinspection PyMethodMayBeStatic
    def web_url(
        self,
        dataset: str,
        namespace: Namespace,
        version: str | None = None,
        sample: str = None,
    ) -> str | None:
        """
        URL corresponding to the WEB bucket.
        """
        url = hail_batch.web_url(
            suffix='', 
            dataset=dataset,
            access_level=namespace.value,
        )
        if version:
            url += f'/{version}'
        if sample:
            url += f'/{sample}'
        return url
