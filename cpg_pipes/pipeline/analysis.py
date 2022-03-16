"""
Sample metadata DB Analysis entry.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import List, Union
from cloudpathlib import CloudPath
from hailtop.batch import ResourceGroup, ResourceFile, Batch

logger = logging.getLogger(__file__)


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
    output: CloudPath | None

    @staticmethod
    def parse(data: dict) -> 'Analysis':
        req_keys = ['id', 'type', 'status']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Analysis" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')
        
        output = data.get('output')
        if output:
            output = CloudPath(output)

        a = Analysis(
            id=int(data['id']),
            type=AnalysisType.parse(data['type']),
            status=AnalysisStatus.parse(data['status']),
            sample_ids=set(data.get('sample_ids', [])),
            output=output,
        )
        return a


class CramPath:
    """
    Represents alignment data on a bucket within the pipeline. 
    Includes a path to a CRAM or a BAM file along with a correponding index,
    and a corresponding fingerprint path.
    """
    def __init__(
        self, 
        path: str | CloudPath, 
        index_path: CloudPath | str | None = None
    ):
        self.path = CloudPath(path)
        self.is_bam = self.path.suffix == '.bam'
        self.ext = 'cram' if not self.is_bam else 'bam'
        self.index_ext = 'crai' if not self.is_bam else 'bai'
        self._index_path = index_path
        self.somalier_path = CloudPath(f'{self.path}.somalier')

    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'CRAM({self.path})'
    
    def exists(self) -> bool:
        return self.path.exists()

    @property
    def index_path(self) -> CloudPath:
        """
        Path to the corresponding index
        """
        return (
            CloudPath(self._index_path) 
            if self._index_path 
            else CloudPath(f'{self.path}.{self.index_ext}')
        )

    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            self.ext: str(self.path),
            f'{self.ext}.{self.index_ext}': str(self.index_path),
        })
    

class GvcfPath:
    """
    Represents GVCF data on a bucket within the pipeline. 
    Includes a path to a GVCF file along with a correponding TBI index,
    and a corresponding fingerprint path.
    """
    def __init__(self, path: CloudPath | str):
        self.path = CloudPath(path)
        self.somalier_path = CloudPath(f'{self.path}.somalier')
        
    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'GVCF({self.path})'

    def exists(self) -> bool:
        return self.path.exists()

    @property
    def tbi_path(self) -> CloudPath:
        """
        Path to the corresponding index
        """
        return CloudPath(f'{self.path}.tbi')
    
    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            'g.vcf.gz': str(self.path),
            'g.vcf.gz.tbi': str(self.tbi_path),
        })


FastqPath = Union[str, CloudPath, ResourceFile]


@dataclass
class FastqPair:
    r1: FastqPath
    r2: FastqPath
    
    def __getitem__(self, i):
        assert i == 0 or i == 1, i
        return [self.r1, self.r2][i]

    def as_resources(self, b) -> 'FastqPair':
        """
        Makes a pair of ResourceFile objects for r1 and r2.
        """
        r1 = self.r1 if isinstance(self.r1, ResourceFile) else b.read_input(str(self.r1))
        r2 = self.r2 if isinstance(self.r2, ResourceFile) else b.read_input(str(self.r2))
        return FastqPair(r1, r2)
    

FastqPairs = List[FastqPair]


# Alignment input can be a CRAM file on a bucket, or a list of Fastq pairs on a bucket.
AlignmentInput = Union[FastqPairs, CramPath]
