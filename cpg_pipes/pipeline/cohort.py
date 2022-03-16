"""
Represents a "cohort" target - all samples from all datasets in the pipeline
"""

import logging

from cpg_pipes.storage import Namespace, StorageProvider
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.target import Target

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Cohort(Target):
    """
    Represents a "cohort" target - all samples from all datasets in the pipeline
    """
    def __init__(self, name: str):
        super().__init__()
        self.name = name
        self._datasets: list[Dataset] = []

    def __repr__(self):
        return self.name

    @property
    def target_id(self) -> str:
        return f'Cohort({self.name}, {len(self.get_datasets())} datasets)'

    def get_datasets(self, only_active: bool = True) -> list[Dataset]:
        """
        Gets list of all datasets. 
        Include only "active" datasets (unless only_active is False)
        """
        return [ds for ds in self._datasets if (ds.active or not only_active)]

    def get_all_samples(self, only_active: bool = True) -> list[Sample]:
        """
        Gets a flat list of all samples from all datasets.
        Include only "active" samples (unless only_active is False)
        """
        all_samples = []
        for proj in self.get_datasets(only_active=only_active):
            all_samples.extend(proj.get_samples(only_active))
        return all_samples

    def get_all_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Gets a flat list of CPG IDs for all samples from all datasets.
        """
        return [s.id for s in self.get_all_samples(only_active=only_active)]

    def add_dataset(
        self, 
        name: str, 
        namespace: Namespace | None = None,
        storage_provider: StorageProvider | None = None,
    ) -> Dataset:
        """
        Create a dataset and add to the cohort.
        """
        ds_by_name = {ds.name: ds for ds in self._datasets}
        if name in ds_by_name:
            logger.warning(f'Dataset {name} already exists in the cohort')
            return ds_by_name[name]
        p = Dataset(
            name=name, 
            namespace=namespace,
            storage_provider=storage_provider,
        )
        self._datasets.append(p)
        return p
