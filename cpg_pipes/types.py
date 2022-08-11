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
from .utils import exists

logger = logging.getLogger(__file__)


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
        index_path: str | Path | None = None,
        reference_assembly: str = None,
    ):
        self.path = to_path(path)
        self.is_bam = self.path.suffix == '.bam'
        self.ext = 'cram' if not self.is_bam else 'bam'
        self.index_ext = 'bai' if self.is_bam else 'crai'
        if not index_path:
            self.full_index_ext = f'{self.ext}.{self.index_ext}'
            self.index_path = self.path.with_suffix(f'.{self.full_index_ext}')
        else:
            assert str(index_path).endswith(self.index_ext)
            self.index_path = to_path(index_path)
            self.full_index_ext = str(self.index_path).replace(self.path.stem, '')
        self.somalier_path = to_path(f'{self.path}.somalier')
        self.reference_assembly = reference_assembly

    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        return f'CRAM({self.path})'

    def exists(self) -> bool:
        """
        CRAM file exists.
        """
        return exists(self.path)

    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(
            **{
                self.ext: str(self.path),
                self.full_index_ext: str(self.index_path),
            }
        )

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
        return FastqPair(
            *[
                self[i]
                if isinstance(self[i], ResourceFile)
                else b.read_input(str(self[i]))
                for i in [0, 1]
            ]
        )

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
        return all(utils.exists(pair.r1) and utils.exists(pair.r2) for pair in self)

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
        return ''.join(
            [
                f'{{{",".join(set(chars))}}}' if len(set(chars)) > 1 else chars[0]
                for chars in zip(*map(str, all_fastq_paths))
            ]
        )
