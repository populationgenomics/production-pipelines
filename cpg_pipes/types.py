"""
Basic bioinformatics status types.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import List, Union

from cpg_pipes import Path, to_path
from hailtop.batch import ResourceGroup, ResourceFile, Batch

logger = logging.getLogger(__file__)


class CramPath:
    """
    Represents alignment data on a bucket within the pipeline. 
    Includes a path to a CRAM or a BAM file along with a correponding index,
    and a corresponding fingerprint path.
    """
    def __init__(
        self, 
        path: str | Path, 
        index_path: Path | str | None = None
    ):
        self.path = to_path(path)
        self.is_bam = self.path.suffix == '.bam'
        self.ext = 'cram' if not self.is_bam else 'bam'
        self.index_ext = 'crai' if not self.is_bam else 'bai'
        self._index_path = index_path
        self.somalier_path = to_path(f'{self.path}.somalier')

    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'CRAM({self.path})'
    
    def exists(self) -> bool:
        """
        CRAM file exists.
        """
        return self.path.exists()

    @property
    def index_path(self) -> Path:
        """
        Path to the corresponding index
        """
        return (
            to_path(self._index_path) 
            if self._index_path 
            else to_path(f'{self.path}.{self.index_ext}')
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
    def __init__(self, path: Path | str):
        self.path = to_path(path)
        self.somalier_path = to_path(f'{self.path}.somalier')
        
    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'GVCF({self.path})'

    def exists(self) -> bool:
        """
        GVCF file exists.
        """
        return self.path.exists()

    @property
    def tbi_path(self) -> Path:
        """
        Path to the corresponding index
        """
        return to_path(f'{self.path}.tbi')
    
    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            'g.vcf.gz': str(self.path),
            'g.vcf.gz.tbi': str(self.tbi_path),
        })


FastqPath = Union[str, Path, ResourceFile]


@dataclass
class FastqPair:
    """
    Pair of FASTQ files
    """
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


class SequencingType(Enum):
    """
    Type (scope) of sequencing experiment.
    """

    WGS = 'wgs'
    EXOME = 'exome'
    SINGLE_CELL = 'single-cell'

    @staticmethod
    def parse(val: str) -> 'SequencingType':
        """
        Parse a string into a SequencingType object.
        """
        d = {v.value: v for v in SequencingType}
        if val not in d:
            raise ValueError(
                f'Unrecognised sequence type {val}. Available: {list(d.keys())}'
            )
        return d[val.lower()]
