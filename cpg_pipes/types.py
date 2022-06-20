"""
Wrappers for bioinformatics file types (CRAM, GVCF, FASTQ, etc).
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Union, Optional

from hailtop.batch import ResourceGroup, ResourceFile, Batch

from . import Path, to_path
from . import utils

logger = logging.getLogger(__file__)


class SequencingType(Enum):
    """
    Type (scope) of a sequencing experiment.
    """

    GENOME = 'genome'
    EXOME = 'exome'
    SINGLE_CELL = 'single_cell'
    MTSEQ = 'mtseq'
    ONT = 'ont'
    
    def seqr_value(self) -> str:
        """
        Map to Seqr-style string
        https://github.com/broadinstitute/seqr/blob/e0c179c36c0f68c892017de5eab2e4c1b9ffdc92/seqr/models.py#L592-L594
        """
        return {
            self.GENOME: 'WGS',
            self.EXOME: 'WES',
            self.SINGLE_CELL: 'RNA',
        }.get(self, '')

    @staticmethod
    def parse(str_val: str) -> 'SequencingType':
        """
        Parse a string into a SequencingType object.
        
        >>> SequencingType.parse('genome')
        SequencingType.GENOME
        >>> SequencingType.parse('wes')
        SequencingType.EXOME
        """
        str_to_val: dict[str, SequencingType] = {} 
        for val, str_vals in {
            SequencingType.GENOME: ['genome', 'wgs'],
            SequencingType.EXOME: ['exome', 'wts', 'wes'],
            SequencingType.SINGLE_CELL: ['single_cell', 'single_cell_rna'],
            SequencingType.MTSEQ: ['mtseq']
        }.items():
            for str_v in str_vals:
                str_v = str_v.lower()
                assert str_v not in str_to_val, (str_v, str_to_val)
                str_to_val[str_v] = val
        str_val = str_val.lower()
        if str_val not in str_to_val:
            raise ValueError(
                f'Unrecognised sequence type {str_val}. '
                f'Available: {list(str_to_val.keys())}'
            )
        return str_to_val[str_val]


class AlignmentInput(ABC):
    """
    Data that works as input for alignment or realignment.
    """
    @abstractmethod
    def exists(self) -> bool:
        """
        Check if all files exist.
        """

    @abstractmethod
    def path_glob(self) -> str:
        """
        Compact representation of file paths. 
        """


class CramPath(AlignmentInput):
    """
    Represents alignment data on a bucket within the pipeline.
    Includes a path to a CRAM or a BAM file along with a corresponding index,
    and a corresponding fingerprint path.
    """

    def __init__(
        self, 
        path: str | Path, 
        index_path: Path | str | None = None,
    ):
        self.path = to_path(path)
        self.is_bam = self.path.suffix == '.bam'
        self.ext = 'cram' if not self.is_bam else 'bam'
        self.index_ext = 'crai' if not self.is_bam else 'bai'
        self._index_path = index_path
        self.somalier_path = to_path(f'{self.path}.somalier')
        self.sequencing_type: Optional[SequencingType] = None

    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'CRAM({self.path})'

    def exists(self) -> bool:
        """
        CRAM file exists.
        """
        return self.path.exists()
    
    def index_exists(self) -> bool:
        """
        CRAI/BAI index exists
        """
        return self.index_path.exists()

    @property
    def index_path(self) -> Path:
        """
        Path to the corresponding CRAI/BAI index
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
        d = {
            self.ext: str(self.path),
        } 
        if self.index_exists():
            d[f'{self.ext}.{self.index_ext}'] = str(self.index_path)

        return b.read_input_group(**d)
    
    def path_glob(self) -> str:
        """
        Compact representation of file paths. 
        For a CRAM, it's just the CRAM file path.
        """
        return str(self.path)


class GvcfPath:
    """
    Represents GVCF data on a bucket within the pipeline.
    Includes a path to a GVCF file along with a corresponding TBI index,
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
        return b.read_input_group(
            **{
                'g.vcf.gz': str(self.path),
                'g.vcf.gz.tbi': str(self.tbi_path),
            }
        )


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
        return FastqPair(*[
            self[i] 
            if isinstance(self[i], ResourceFile) 
            else b.read_input(str(self[i]))
            for i in [0, 1]
        ])

    def __str__(self):
        return f'{self.r1}|{self.r2}'


class FastqPairs(list[FastqPair], AlignmentInput):
    """
    Multiple FASTQ file pairs belonging to the same sample 
    (e.g. multiple lanes or top-ups).
    """
    def exists(self) -> bool:
        """
        Check if each FASTQ file in each pair exist.
        """
        return all(
            utils.exists(pair.r1) and utils.exists(pair.r2) 
            for pair in self
        )

    def path_glob(self) -> str:
        """
        Compact representation of file paths. 
        For FASTQ pairs, it's glob string to find all FASTQ files.

        >>> FastqPairs([
        >>>     FastqPair('gs://sample_R1.fq.gz', 'gs://sample_R2.fq.gz'),
        >>> ]).path_glob()
        'gs://sample_R{2,1}.fq.gz'
        >>> FastqPairs([
        >>>     FastqPair('gs://sample_L1_R1.fq.gz', 'gs://sample_L1_R2.fq.gz'), 
        >>>     FastqPair('gs://sample_L2_R1.fq.gz', 'gs://sample_L2_R2.fq.gz'),
        >>> ]).path_glob()
        'gs://sample_L{2,1}_R{2,1}.fq.gz'
        """
        all_fastq_paths = []
        for pair in self:
            all_fastq_paths.extend([pair.r1, pair.r2])
        # Triple braces are intentional: they are resolved to single ones.
        return ''.join([
            f'{{{",".join(set(chars))}}}' if len(set(chars)) > 1 else chars[0] 
            for chars in zip(*map(str, all_fastq_paths))
        ])
