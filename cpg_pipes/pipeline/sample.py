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


@dataclass(init=False)
class Sample(Target):
    """
    Corresponds to one Sample entry in the SMDB
    """
    def __init__(
        self, 
        id: str, 
        external_id: str, 
        project,
        participant_id: Optional[str] = None,
        meta: dict = None,
    ):
        super().__init__()
        self.id = id
        self.external_id = external_id
        self.project = project
        self.participant_id = participant_id or external_id
        self.meta = meta or dict()
        self.alignment_input: Optional[AlignmentInput] = None
        self.seq: Optional[SmSequence] = None
        self.pedigree: Optional[PedigreeInfo] = None

    @property
    def unique_id(self) -> str:
        return self.id

    def get_ped_dict(self, use_participant_id: bool = False) -> Dict:
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
