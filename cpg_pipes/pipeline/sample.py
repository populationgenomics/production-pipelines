"""
Corresponds to one Sample entry in the SMDB
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from cpg_pipes.pipeline.analysis import CramPath, GvcfPath, AlignmentInput
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sequence import SmSequence

logger = logging.getLogger(__file__)


class Sample(Target):
    """
    Represents a Sample
    """
    def __init__(
        self, 
        id: str,  # pylint: disable=redefined-builtin
        external_id: str,
        dataset: 'Dataset',  # type: ignore  # noqa: F821
        participant_id: str | None = None,
        meta: dict | None = None,
        seq: SmSequence | None = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input: AlignmentInput | None = None
    ):
        super().__init__()
        self.id = id
        self.external_id = external_id
        self.dataset = dataset
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.seq = seq
        self.pedigree: PedigreeInfo | None
        if pedigree:
            self.pedigree = pedigree
        elif 'sex' in self.meta:
            self.pedigree = PedigreeInfo(
                sample=self,
                fam_id=self.participant_id,
                sex=Sex.parse(self.meta['sex'])
            )
        else:
            self.pedigree = None
        self._alignment_input = alignment_input

    def __repr__(self):
        return (
            f'Sample({self.dataset.name}/{self.id}|{self.external_id}' +
            (f', participant={self._participant_id}' 
             if self._participant_id else '') +
            f', forced={self.forced}' +
            f', active={self.active}' +
            f', meta={self.meta}' +
            (f', seq={self.seq}' if self.seq else '') +
            (f', alignment_input={self._alignment_input}' 
             if self._alignment_input else '') +
            (f', pedigree={self.pedigree}' if self.pedigree else '') +
            f')'
        )

    def __str__(self):
        ai_tag = ''
        if self._alignment_input:
            if isinstance(self._alignment_input, CramPath):
                if self._alignment_input.is_bam:
                    ai_tag = f'|SEQ=CRAM'
                else:
                    ai_tag = f'|SEQ=BAM'
            else:
                ai_tag = f'|SEQ={len(self._alignment_input)}FQS'

        return f'Sample({self.dataset.name}/{self.id}|{self.external_id}{ai_tag})'

    @property
    def alignment_input(self) -> AlignmentInput | None:
        """
        Returns input for (re-)alignment.
        """
        if self._alignment_input:
            return self._alignment_input

        if self.seq:
            return self.seq.alignment_input
        return None

    @alignment_input.setter
    def alignment_input(self, value: AlignmentInput):
        """
        Set alignment_input
        """
        self._alignment_input = value

    @property
    def participant_id(self) -> str:
        """
        Get participant's ID. 
        Uses external_id whenever participant_id is not available in the DB
        """
        return self._participant_id or self.external_id

    @property
    def target_id(self) -> str:
        return f'Sample({self.id})'

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
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

    def get_cram_path(self) -> CramPath:
        """
        Path to a CRAM file. Not checking its existence here.
        """
        return CramPath(self.dataset.get_bucket() / 'cram' / f'{self.id}.cram')

    def get_gvcf_path(self) -> GvcfPath:
        """
        Path to a GVCF file. Not checking its existence here.
        """
        return GvcfPath(self.dataset.get_bucket() / 'gvcf' / f'{self.id}.g.vcf.gz')


class Sex(Enum):
    """
    Sex as in PED format
    """
    UNKNOWN = 0
    MALE = 1
    FEMALE = 2
    
    @staticmethod
    def parse(sex: str) -> 'Sex':
        if sex.lower() in ('m', 'male', '1'):
            return Sex.MALE
        if sex.lower() in ('f', 'female', '2'):
            return Sex.FEMALE
        return Sex.UNKNOWN


@dataclass
class PedigreeInfo:
    """
    Pedigree relationsips with other samples in the cohort, and other PED data
    """
    sample: Sample
    sex: Sex
    fam_id: str | None = None
    phenotype: str | None = None
    dad: Sample | None = None
    mom: Sample | None = None

    def get_ped_dict(self, use_participant_id: bool = False) -> dict:
        """
        Returns a dictionary of pedigree fields for this sample, corresponging
        a PED file entry.
        """
        def _get_id(_s: Sample|None):
            if _s is None:
                return '0'
            if use_participant_id:
                return _s.participant_id
            return _s.id

        return {
            'Family.ID': self.fam_id or self.sample.participant_id,
            'Individual.ID': _get_id(self.sample),
            'Father.ID': _get_id(self.dad),
            'Mother.ID': _get_id(self.mom),
            'Sex': str(self.sex.value),
            'Phenotype': self.phenotype,
        }
