"""
Sample metadata DB Analysis entry.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import cast

from hailtop.batch import ResourceGroup, ResourceFile, Batch

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

    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            self.ext: self.path,
            f'{self.ext}.{self.index_ext}': self.index_path,
        })
    
    def alignment_input(self) -> 'AlignmentInput':
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
    
    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            'g.vcf.gz': self.path,
            'g.vcf.gz.tbi': self.tbi_path,
        })


class AlignmentInput:
    """
    Represents inputs for an alignment job, which can be a set of fastq files,
    or a CRAM or a BAM file with an index.
    """
    def __init__(
        self, 
        fqs1: list[str|ResourceFile]|None = None,
        fqs2: list[str|ResourceFile]|None = None,
        cram_path: CramPath|Path|str|ResourceGroup|None = None,
        index_path: Path|str|None = None,
    ):
        self.fqs1 = fqs1
        self.fqs2 = fqs2
        self.cram_path: CramPath|ResourceGroup|None = None
        if isinstance(cram_path, str|Path):
            self.cram_path = CramPath(cram_path, index_path=index_path)
        else:
            self.cram_path = cram_path

    def __repr__(self):
        return (
            f'AlignmentInput(' +
            (self.cram_path if self.is_bam_or_cram() else str((self.fqs1, self.fqs2))) +
            f')'
        )

    def is_fastq(self) -> bool:
        """
        Checks that it's a fastq pair, and both in pair are of the same type and length
        """
        if self.fqs1 or self.fqs2:
            assert self.fqs1 and self.fqs2, self
            if any(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]), self
            elif any(isinstance(fq, ResourceFile) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, ResourceFile) for fq in [self.fqs1, self.fqs2]), self
            else:
                assert len(self.fqs1) == len(self.fqs2), self
            return True
        assert self.cram_path, self
        return False

    def is_bam_or_cram(self) -> bool:
        """
        Checks that it's a BAM or a CRAM file
        """
        if self.cram_path:
            return True
        assert self.fqs1 and self.fqs2, self
        return False

    def get_fqs1(self) -> list[str|ResourceFile]:
        assert self.is_fastq()
        return cast(list, self.fqs1)

    def get_fqs2(self) -> list[str|ResourceFile]:
        assert self.is_fastq()
        return cast(list, self.fqs2)

    def as_fq_inputs(self, b) -> tuple[list[ResourceFile], list[ResourceFile]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        assert self.is_fastq()
        self.fqs1 = cast(list, self.fqs1)
        self.fqs2 = cast(list, self.fqs2)
        if isinstance(self.fqs1[0], ResourceFile):
            files1 = self.fqs1
            files2 = self.fqs2
        else:
            files1 = [b.read_input(f1) for f1 in self.fqs1]
            files2 = [b.read_input(f1) for f1 in self.fqs2]
        return files1, files2

    def as_cram_input_group(self, b) -> ResourceGroup:
        """
        Makes a ResourceGroup of bam/cram with accompanying index
        """
        assert self.is_bam_or_cram()

        if isinstance(self.cram_path, ResourceGroup):
            return cast(ResourceGroup, self.cram_path)

        self.cram_path = cast(CramPath, self.cram_path)
        return self.cram_path.resource_group(b)
