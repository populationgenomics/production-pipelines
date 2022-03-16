"""
Represents a CPG dataset in a particular namespace: main and test.
"""

import logging
from cloudpathlib import CloudPath

from cpg_pipes.storage import Namespace, StorageProvider
from cpg_pipes.pipeline.sequence import SmSequence
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sample import Sample, PedigreeInfo

logger = logging.getLogger(__file__)


class Dataset(Target):
    """
    Represents a CPG dataset in a particular namespace: main or test.

    Each `dataset` at the CPG corresponds to
    * one GCP project: https://github.com/populationgenomics/team-docs/tree/main/storage_policies
    * one Pulumi stack: https://github.com/populationgenomics/analysis-runner/tree/main/stack
    * two sample metadata projects: main and test (the latter has a `-test` ending).
    
    An object of this class is parametrised by a dataset name and a namespace,
    meaning that it matches exactly one GCP project, exactly one stack, and exactly 
    one sample metadata project.

    An object has two ID-like fields: `stack` and `name`:
    * `stack` is the name of the dataset (matches names of a GCP project or
       a Pulumi stack), e.g. "seqr", "hgdp".
    * `name` is the name of the namespace-specific sample-metadata project,
       e.g. "seqr", "seqr-test", "hgdp", "hgdp-test".
    """
    def __init__(
        self, 
        name: str,
        namespace: Namespace | None = None,
        storage_provider: StorageProvider | None = None,
    ):
        """
        Input `name` can be either e.g. "seqr" or "seqr-test". The latter will be 
        resolved to stack="seqr" and is_test=True, unless `namespace` is provided 
        explicitly.

        Also note that a Pipeline is passed as a parameter. A Dataset can exist 
        outside of a cohort, e.g. if it's an analysis dataset that exists only 
        to track joint analysis of multiple datasets, but itself doesn't contain 
        any samples.
        """
        super().__init__()
        
        self._storage_provider = storage_provider

        self._samples: list[Sample] = []

        self.namespace = namespace or Namespace.MAIN
        if name.endswith('-test'):
            self.stack = name[:-len('-test')]
            if self.namespace == Namespace.MAIN:
                self.namespace = Namespace.TEST
        else:
            self.stack = name

    @property
    def is_test(self) -> bool:
        return self.namespace != Namespace.MAIN
    
    @property
    def name(self) -> str:
        return self.stack + ('-test' if self.is_test else '')

    @property
    def target_id(self) -> str:
        return f'Dataset({self.name})'
    
    def __repr__(self):
        return f'Dataset("{self.name}", {len(self.get_samples())} samples)'

    def __str__(self):
        return f'{self.name} ({len(self.get_samples())} samples)'

    @property
    def storage_provider(self) -> StorageProvider:
        if not self._storage_provider:
            raise ValueError(
                '_storage_provider is not set. storage_provider must be passed to the '
                'Dataset() constructor before calling Dataset.get_bucket()'
            )
        return self._storage_provider

    def get_bucket(self, **kwargs) -> CloudPath:
        """
        The primary dataset bucket (-main or -test).
        """
        return self.storage_provider.get_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs,
        )

    def get_tmp_bucket(self, **kwargs) -> CloudPath:
        """
        The tmp bucket (-main-tmp or -test-tmp)
        """
        return self.storage_provider.get_tmp_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def get_analysis_bucket(self, **kwargs) -> CloudPath:
        """
        Get analysis bucket (-main-analysis or -test-analysis)
        """
        return self.storage_provider.get_analysis_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def get_web_bucket(self, **kwargs) -> CloudPath:
        """
        Get web bucket (-main-web or -test-web)
        """
        return self.storage_provider.get_web_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )

    def get_web_url(self, **kwargs) -> CloudPath:
        """
        Get web base URL
        """
        return self.storage_provider.get_web_url(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def add_sample(
        self, 
        id: str,  # pylint: disable=redefined-builtin
        external_id: str, 
        participant_id: str | None = None,
        seq: SmSequence | None = None,
        pedigree: PedigreeInfo | None = None,
        **kwargs
    ) -> Sample:
        """
        Create a new sample and add it into the dataset.
        """
        s = Sample(
            id=id, 
            external_id=external_id,
            participant_id=participant_id,
            seq=seq,
            pedigree=pedigree,
            dataset=self,
            meta=kwargs,
        )
        self._samples.append(s)
        return s
    
    def get_samples(self, only_active: bool = True) -> list[Sample]:
        """
        Get dataset's samples. Inlcude only "active" samples, unless only_active=False
        """
        return [s for s in self._samples if (s.active or not only_active)]

    def get_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Get dataset's sample IDs. Inlcude only "active" samples, 
        unless only_active=False.
        """
        return [s.id for s in self.get_samples(only_active=only_active)]
