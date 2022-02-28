"""
Corresponds to one Sample entry in the SMDB
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional

from cpg_pipes.hb.inputs import AlignmentInput
from cpg_pipes.pipeline.target import Target
from cpg_pipes.smdb.types import SmSequence

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Sample(Target):
    """
    Represents a Sample
    """
    def __init__(
        self, 
        id: str,  # pylint: disable=redefined-builtin
        external_id: str,
        project: 'Project',  # type: ignore  # noqa: F821
        participant_id: str | None = None,
        meta: dict | None = None
    ):
        super().__init__()
        self.id = id
        self.external_id = external_id
        self.project = project
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.alignment_input: Optional[AlignmentInput] = None
        self.seq: Optional[SmSequence] = None
        self.pedigree: Optional['PedigreeInfo'] = None

    def __repr__(self):
        return (
            f'Sample({self.id}|{self.external_id}' +
            (f'participant_id={self._participant_id}, ' 
             if self._participant_id else '') +
            f', project={self.project.id}'
            f', forced={self.forced}'
            f', active={self.active}'
            f', alignment_input={self.alignment_input}'
            f', meta={self.meta}'
            f', seq={self.seq}'
            f', pedigree={self.pedigree}'
            f')'
        )

    @property
    def participant_id(self) -> str:
        """
        Get participant's ID. 
        Uses external_id whenever participant_id is not available in the DB
        """
        return self._participant_id or self.external_id

    @property
    def unique_id(self) -> str:
        return self.id

    def get_ped_dict(self, use_participant_id: bool = False) -> Dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sample, corresponging
        a PED file entry.
        """
        if self.pedigree:
            return self.pedigree.get_ped_dict(use_participant_id)
        return {
            'Family.ID': self.participant_id if use_participant_id else self.id,
            'Individual.ID': self.participant_id if use_participant_id else self.id,
            'Father.ID': '0',
            'Mother.ID': '0',
            'Sex': '0',
            'Phenotype': '0',
        }


class Sex(Enum):
    """
    Sex as in PED format
    """
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2


@dataclass
class PedigreeInfo:
    """
    Pedigree relationsips with other samples in the cohort, and other PED data
    """
    sample: Sample
    fam_id: str
    dad: Optional[Sample]
    mom: Optional[Sample]
    sex: Sex
    phenotype: str

    def get_ped_dict(self, use_participant_id: bool = False) -> Dict:
        """
        Returns a dictionary of pedigree fields for this sample, corresponging
        a PED file entry.
        """
        def _get_id(_s: Optional[Sample]):
            if _s is None:
                return '0'
            if use_participant_id:
                return _s.participant_id
            return _s.id

        return {
            'Family.ID': self.fam_id,
            'Individual.ID': _get_id(self.sample),
            'Father.ID': _get_id(self.dad),
            'Mother.ID': _get_id(self.mom),
            'Sex': str(self.sex.value),
            'Phenotype': self.phenotype,
        }
