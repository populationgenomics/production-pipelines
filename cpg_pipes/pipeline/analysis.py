"""
Sample metadata DB Analysis entry.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import hailtop.batch as hb
from hailtop.batch import ResourceGroup

from cpg_pipes.alignment_input import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class AnalysisType(Enum):
    """
    Corresponds to SMDB Analysis types:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L4-L11

    Re-defined in a separate module to decouple from the main sample-metadata module,
    so decorators can use `@stage(analysis_type=AnalysisType.QC)` without importing
    the sample-metadata package.
    """
    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'
    
    @staticmethod
    def parse(name: str) -> 'AnalysisType':
        return {v.value: v for v in AnalysisType}[name.lower()]


class AnalysisStatus(Enum):
    """
    Corresponds to SMDB Analysis statuses:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L14-L21
    """
    QUEUED = 'queued'
    IN_PROGRESS = 'in-progress'
    FAILED = 'failed'
    COMPLETED = 'completed'
    UNKNOWN = 'unknown'
    
    @staticmethod
    def parse(name: str) -> 'AnalysisStatus':
        return {v.value: v for v in AnalysisStatus}[name.lower()]


@dataclass
class Analysis:
    """
    Sample metadata DB Analysis entry.

    See the sample-metadata package for more details: 
    https://github.com/populationgenomics/sample-metadata
    """
    id: int
    type: AnalysisType
    status: AnalysisStatus
    sample_ids: set[str]
    output: Path|None

    @staticmethod
    def parse(data: dict) -> 'Analysis':
        req_keys = ['id', 'type', 'status']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Analysis" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')
        
        type_ = AnalysisType.parse(data['type'])
        
        output = data.get('output')
        if output:
            output = Path(output)

        a = Analysis(
            id=int(data['id']),
            type=type_,
            status=data['status'],
            sample_ids=set(data.get('sample_ids', [])),
            output=output,
        )
        return a


class CramPath:
    """
    Represents a path to CRAM or a BAM file along with a correponding index,
    and a corresponding fingerprint path.
    """
    def __init__(self, path: str|Path, index_path: Path|str|None = None):
        self.path: Path = Path(path)
        self.is_bam = self.path.suffix == '.bam'
        self.ext = 'cram' if not self.is_bam else 'bam'
        self.index_ext = 'crai' if not self.is_bam else 'bai'
        self._index_path = index_path
        self.somalier_path = Path(self.path.stem + '.somalier')

    @property
    def index_path(self) -> Path:
        """
        Path to the corresponding index
        """
        return (
            Path(self._index_path) 
            if self._index_path 
            else Path(f'{self.path}.{self.index_ext}')
        )

    def resource_group(self, b: hb.Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            self.ext: self.path,
            f'{self.ext}.{self.index_ext}': self.index_path,
        })
    
    def alignment_input(self) -> AlignmentInput:
        return AlignmentInput(cram_path=self)


class GvcfPath:
    """
    Represents a GVCF file path alon g with a corresponding index, 
    and a corresponding fingerprint path.
    """
    def __init__(self, path: Path|str):
        self.path = Path(path)
        self.somalier_path = Path(self.path.stem + '.somalier')

    @property
    def tbi_path(self) -> Path:
        """
        Path to the corresponding index
        """
        return Path(f'{self.path}.tbi')
    
    def resource_group(self, b: hb.Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            'g.vcf.gz': self.path,
            'g.vcf.gz.tbi': self.tbi_path,
        })
