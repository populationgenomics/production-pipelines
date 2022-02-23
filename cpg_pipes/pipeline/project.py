"""
A CPG stack (aka dataset), or a project from the sample-metadata database.
"""
from typing import List, Optional
import logging

from cpg_pipes.namespace import Namespace
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sample import Sample

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Project(Target):
    """
    Represents a CPG stack (aka dataset), or a project from sample-metadata database.
    """
    def __init__(
        self, 
        name: str,
        pipeline,
        namespace: Optional[Namespace] = None,
    ):
        """
        Can have a "name" and "stack":
        * "name" is in the SMDB terms: e.g. can be "seqr", "seqr-test".
        * "stack" is in the CPG storage policy terms, so can be only e.g. "seqr", but not "seqr-test"
        """
        super().__init__()
        if name.endswith('-test'):
            self.is_test = True
            self.stack = name[:-len('-test')]
        else:
            self.is_test = False
            self.stack = name
            
        self.name = self.stack
        if self.is_test:
            self.name = self.stack + '-test'

        self.pipeline = pipeline
        self.is_test = namespace != Namespace.MAIN
        self._samples: List[Sample] = []

    def get_bucket(self):
        """
        The primary project bucket (-main or -test) 
        """
        return f'gs://cpg-{self.stack}-{self.pipeline.output_suf}'

    def get_tmp_bucket(self):
        """
        The tmp bucket (-main-tmp or -test-tmp)
        """
        return (
            f'gs://cpg-{self.stack}-{self.pipeline.output_suf}-tmp/'
            f'{self.pipeline.name}/'
            f'{self.pipeline.output_version}'
        )

    def __repr__(self):
        return self.name

    @property
    def unique_id(self) -> str:
        return self.name

    def add_sample(
        self, 
        id: str, 
        external_id: str, 
        participant_id: Optional[str] = None,
        **kwargs
    ) -> Sample:
        s = Sample(
            id=id, 
            external_id=external_id,
            participant_id=participant_id,
            project=self,
            meta=kwargs,
        )
        self._samples.append(s)
        return s
    
    def get_samples(self, only_active: bool = True) -> List[Sample]:
        """
        Get project's samples. Inlcude only "active" samples, unless only_active
        """
        return [s for s in self._samples if (s.active or not only_active)]
