import logging
from dataclasses import dataclass, field, InitVar
from enum import Enum
from typing import Dict, Optional

from cpg_pipes.hb.inputs import AlignmentInput
from cpg_pipes.pipeline.target import Target
from cpg_pipes.smdb.types import SmSequence

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@dataclass
class Sample(Target):
    """
    Corresponds to one Sample entry in the SMDB
    """
    id: str
    external_id: str
    project: 'Project'  # type: ignore  # noqa: F821
    participant_id: InitVar[Optional[str]] = None
    meta: dict = field(default_factory=dict)
    alignment_input: Optional[AlignmentInput] = None
    seq: Optional[SmSequence] = field(repr=False, default=None)
    pedigree: Optional['PedigreeInfo'] = None
    
    def __post_init__(self, participant_id: Optional[str]):
        self.participant_id: str = participant_id or self.external_id

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
        else:
            return {
                'Family.ID': self.participant_id if use_participant_id else self.id,
                'Individual.ID': self.participant_id if use_participant_id else self.id,
                'Father.ID': '0',
                'Mother.ID': '0',
                'Sex': '0',
                'Phenotype': '0',
            }


class Sex(Enum):
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2


@dataclass
class PedigreeInfo:
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
            elif use_participant_id:
                return _s.participant_id
            else:
                return _s.id

        return {
            'Family.ID': self.fam_id,
            'Individual.ID': _get_id(self.sample),
            'Father.ID': _get_id(self.dad),
            'Mother.ID': _get_id(self.mom),
            'Sex': str(self.sex.value),
            'Phenotype': self.phenotype,
        }
