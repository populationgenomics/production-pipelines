"""
Wrappers for bioinformatics file types (CRAM, GVCF, FASTQ, etc).
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Union

from hailtop.batch import Batch, ResourceFile, ResourceGroup

from cpg_utils import Path, to_path

from .utils import exists


class AlignmentInput(ABC):
    """
    Data that works as input for alignment or realignment.
    """

    @abstractmethod
    def exists(self) -> bool:
        """
        Check if all files exist.
        """


class CramOrBamPath(AlignmentInput, ABC):
    """
    Represents a path to a CRAM or a BAM file, optionally with corresponding index.
    """

    def __init__(
        self,
        path: str | Path,
        index_path: str | Path | None = None,
        reference_assembly: str | Path | None = None,
    ):
        self.path = to_path(path)
        self.index_path: Path | None = None
        self.full_index_suffix: str | None = None
        if index_path:
            self.index_path = to_path(index_path)
            assert self.index_path.suffix == f'.{self.index_ext}'
            self.full_index_suffix = str(self.index_path).replace(str(self.path.with_suffix('')), '')
        self.reference_assembly = None
        if reference_assembly:
            self.reference_assembly = to_path(reference_assembly)

    @property
    @abstractmethod
    def ext(self) -> str:
        """The canonical extension for the file type, without a '.' at the start."""

    @property
    @abstractmethod
    def index_ext(self) -> str:
        """The canonical index file extension, without a '.' at the start."""

    def __str__(self) -> str:
        return str(self.path)

    def __repr__(self) -> str:
        """
        >>> repr(CramPath('gs://bucket/sequencing_group.cram', 'gs://bucket/sequencing_group.cram.crai'))
        'CRAM(gs://bucket/sequencing_group.cram+.cram.crai)'
        """
        res = str(self.path)
        if self.index_path:
            assert self.full_index_suffix
            res += f'+{self.full_index_suffix}'
        return f'{self.ext.upper()}({res})'

    def exists(self) -> bool:
        """
        CRAM file exists.
        """
        return exists(self.path)

    def resource_group(self, b: Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        d = {
            self.ext: str(self.path),
        }
        if self.full_index_suffix:
            d[self.full_index_suffix] = str(self.index_path)

        return b.read_input_group(**d)


class BamPath(CramOrBamPath):
    """
    Represents a path to a BAM file, optionally with corresponding index.
    """

    EXT = 'bam'
    INDEX_EXT = 'bai'

    def __init__(
        self,
        path: str | Path,
        index_path: str | Path | None = None,
    ):
        super().__init__(path, index_path)

    @property
    def ext(self) -> str:
        return BamPath.EXT

    @property
    def index_ext(self) -> str:
        return BamPath.INDEX_EXT


class CramPath(CramOrBamPath):
    """
    Represents a path to a CRAM file, optionally with corresponding index.
    """

    EXT = 'cram'
    INDEX_EXT = 'crai'

    def __init__(
        self,
        path: str | Path,
        index_path: str | Path | None = None,
        reference_assembly: str | Path | None = None,
    ):
        super().__init__(path, index_path, reference_assembly)
        self.somalier_path = to_path(f'{self.path}.somalier')

    @property
    def ext(self) -> str:
        return CramPath.EXT

    @property
    def index_ext(self) -> str:
        return CramPath.INDEX_EXT


class GvcfPath:
    """
    Represents GVCF data on a bucket within the workflow.
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
            },
        )


FastqPath = Union[str, Path, ResourceFile]


@dataclass
class FastqPair(AlignmentInput):
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
            *[(self[i] if isinstance(self[i], ResourceFile) else b.read_input(str(self[i]))) for i in [0, 1]],
        )

    def exists(self) -> bool:
        """
        Check if each FASTQ file in the pair exists.
        """
        return exists(self.r1) and exists(self.r2)

    def __repr__(self):
        """
        Glob string to find all FASTQ files.

        >>> str(FastqPair('gs://sequencing_group_R1.fq.gz', 'gs://sequencing_group_R2.fq.gz'))
        'gs://sequencing_group_R{1,2}.fq.gz'
        """
        return ''.join(
            f'{{{",".join(sorted(set(chars)))}}}' if len(set(chars)) > 1 else chars[0]
            for chars in zip(str(self.r1), str(self.r2))
        )


class FastqPairs(list[FastqPair], AlignmentInput):
    """
    Multiple FASTQ file pairs belonging to the same sequencing_group
    (e.g. multiple lanes or top-ups).
    """

    def exists(self) -> bool:
        """
        Check if each FASTQ file in each pair exists.
        """
        return all(pair.exists() for pair in self)

    def __repr__(self) -> str:
        """
        Glob string to find all FASTQ files.

        >>> repr(FastqPairs([FastqPair('gs://sequencing_group_R1.fq.gz', 'gs://sequencing_group_R2.fq.gz')]))
        'gs://sequencing_group_R{1,2}.fq.gz'
        >>> p1 = FastqPair('gs://sequencing_group_L2_R1.fq.gz', 'gs://sequencing_group_L2_R2.fq.gz')
        >>> p2 = FastqPair('gs://sequencing_group_L1_R1.fq.gz', 'gs://sequencing_group_L1_R2.fq.gz')
        >>> repr(FastqPairs([p1, p2]))
        'gs://sequencing_group_L{1,2}_R{1,2}.fq.gz'
        """
        return ''.join(
            f'{{{",".join(sorted(set(chars)))}}}' if len(set(chars)) > 1 else chars[0]
            for chars in zip(*[repr(pair) for pair in self])
        )
