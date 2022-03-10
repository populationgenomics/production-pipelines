"""
Represents a CPG dataset in a particular namespace: main and test.
"""

from typing import List
import logging

from cloudpathlib import CloudPath

from cpg_pipes.storage import Namespace
from cpg_pipes.pipeline.sequence import SmSequence
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sample import Sample, PedigreeInfo

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


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
        pipeline: 'Pipeline',  # type: ignore  # noqa: F821
        namespace: Namespace | None = None,
    ):
        """
        Input `name` can be either e.g. "seqr" or "seqr-test". The latter will be 
        resolved to stack="seqr" and is_test=True, unless `namespace` is provided 
        explicitly.

        Also note that a Pipeline is passed as a parameter. A Dataset can exist outside
        of a Cohort, e.g. if it's an analysis dataset that exists only to track
        joint analysis of multiple datasets, but itself doesn't contain any samples.
        """
        super().__init__()

        self.pipeline = pipeline
        self._samples: List[Sample] = []

        if name.endswith('-test'):
            self.is_test = True
            self.stack = name[:-len('-test')]
        else:
            self.is_test = False
            self.stack = name
            
        if namespace is not None:
            self.is_test = namespace != Namespace.MAIN

    @property
    def name(self) -> str:
        return self.stack + ('-test' if self.is_test else '')

    @property
    def unique_id(self) -> str:
        return self.name

    def __repr__(self):
        return f'Dataset("{self.name}", {len(self.get_samples())} samples)'

    def __str__(self):
        return f'{self.name} ({len(self.get_samples())} samples)'

    def get_bucket(self) -> CloudPath:
        """
        The primary dataset bucket (-main or -test) 
        """
        prefix = self.pipeline.storage_provider.value
        return (
            CloudPath(f'{prefix}://cpg-{self.stack}-{self.pipeline.output_suf}')
        )

    def get_tmp_bucket(self) -> CloudPath:
        """
        The tmp bucket (-main-tmp or -test-tmp)
        """
        prefix = self.pipeline.storage_provider.value
        return (
            CloudPath(f'{prefix}://cpg-{self.stack}-{self.pipeline.output_suf}-tmp')
            / self.pipeline.name 
            / self.pipeline.output_version
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
