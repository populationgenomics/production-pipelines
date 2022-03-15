"""
Corresponds to one Sample entry in the SMDB
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from cpg_pipes.pipeline.analysis import CramPath, GvcfPath, AnalysisType, AlignmentInput
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sequence import SmSequence

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
            f'Sample({self.id}|{self.external_id}' +
            (f', participant_id={self._participant_id}' 
             if self._participant_id else '') +
            f', dataset={self.dataset.name}' +
            f', forced={self.forced}' +
            f', active={self.active}' +
            f', meta={self.meta}' +
            (f', seq={self.seq}' if self.seq else '') +
            (f', alignment_input={self._alignment_input}' 
             if self._alignment_input else '') +
            (f', pedigree={self.pedigree}' if self.pedigree else '') +
            f')'
        )

    @property
    def alignment_input(self) -> AlignmentInput | None:
        """
        Returns input for (re-)alignment.
        :param use_cram_analysis: check CRAM analysis entry.
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
    def unique_id(self) -> str:
        return self.id

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

    def get_cram_path(self, pipeline) -> CramPath:
        """
        Path to a CRAM file. Not checking its existence here.
        """
        return CramPath(f'{self.dataset.get_bucket(pipeline)}/cram/{self.id}.cram')

    def get_gvcf_path(self, pipeline) -> GvcfPath:
        """
        Path to a GVCF file. Not checking its existence here.
        """
        return GvcfPath(f'{self.dataset.get_bucket(pipeline)}/gvcf/{self.id}.g.vcf.gz')

    def analysis_cram_path(self) -> CramPath | None:
        """
        Get CramPath from analysis SMDB entry.
        """
        if analysis := self.analysis_by_type.get(AnalysisType.CRAM):
            return CramPath(analysis.output)
        return None


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
