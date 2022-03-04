"""
Converts paths into Hail Batch inputs: ResourceFile and ResourceGroup objects.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union, List, Tuple, cast

import hailtop.batch as hb
from hailtop.batch import ResourceGroup

from cpg_pipes import ref_data


def fasta_group(b: hb.Batch) -> ResourceGroup:
    """
    Returns fasta reference resource group
    """
    return b.read_input_group(**ref_data.REF_D)


class Cram:
    """
    Represents a CRAM or a BAM file with a correponding index,
    and a corresponding fingerprint path.
    """
    def __init__(self, path: Path|str, index_path: Path|None = None):
        self.path = Path(path)
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
        return self._index_path or Path(f'{self.path}.{self.index_ext}')

    def __repr__(self) -> str:
        name = 'CRAM' if not self.is_bam else 'BAM'
        return f'{name}({self.path})'
    
    def resource_group(self, b: hb.Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            self.ext: self.path,
            f'{self.ext}.{self.index_ext}': self.index_path,
        })


class Gvcf:
    """
    Represents a GVCF file with a corresponding index, and a corresponding 
    fingerprint path.
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

    def __repr__(self) -> str:
        return f'GVCF({self.path})'
    
    def resource_group(self, b: hb.Batch) -> ResourceGroup:
        """
        Create a Hail Batch resource group
        """
        return b.read_input_group(**{
            'g.vcf.gz': self.path,
            'g.vcf.gz.tbi': self.tbi_path,
        })


@dataclass
class AlignmentInput:
    """
    Represents inputs for an alignment job, which can be a set of fastq files,
    or a CRAM or a BAM file with an index.
    """
    fqs1: Optional[List[Union[str, hb.ResourceFile]]] = None
    fqs2: Optional[List[Union[str, hb.ResourceFile]]] = None
    bam_or_cram_path: Optional[Union[str, hb.ResourceGroup]] = None
    index_path: Optional[str] = None
    
    def is_fastq(self) -> bool:
        """
        Checks that it's a fastq pair, and both in pair are of the same type and length
        """
        if self.fqs1 or self.fqs2:
            assert self.fqs1 and self.fqs2, self
            if any(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]), self
            elif any(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]), self
            else:
                assert len(self.fqs1) == len(self.fqs2), self
            return True
        assert self.bam_or_cram_path, self
        return False

    def is_bam_or_cram(self) -> bool:
        """
        Checks that it's a BAM or a CRAM file
        """
        if self.bam_or_cram_path:
            return True
        assert self.fqs1 and self.fqs2, self
        return False

    def get_fqs1(self) -> List[Union[str, hb.ResourceFile]]:
        assert self.is_fastq()
        return cast(List, self.fqs1)

    def get_fqs2(self) -> List[Union[str, hb.ResourceFile]]:
        assert self.is_fastq()
        return cast(List, self.fqs2)

    def as_fq_inputs(self, b) -> Tuple[List[hb.Resource], List[hb.Resource]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        assert self.is_fastq()
        self.fqs1 = cast(List, self.fqs1)
        self.fqs2 = cast(List, self.fqs2)
        if isinstance(self.fqs1[0], hb.Resource):
            files1 = self.fqs1
            files2 = self.fqs2
        else:
            files1 = [b.read_input(f1) for f1 in self.fqs1]
            files2 = [b.read_input(f1) for f1 in self.fqs2]
        return files1, files2

    def as_cram_input_group(self, b) -> hb.ResourceGroup:
        """
        Makes a ResourceGroup of bam/cram with accompanying index
        """
        assert self.is_bam_or_cram()
        self.bam_or_cram_path = cast(str, self.bam_or_cram_path)
        index_path = self.index_path
        if not index_path:
            if self.bam_or_cram_path.endswith('.bam'):
                index_path = self.bam_or_cram_path + '.bai'
            else:
                assert self.bam_or_cram_path.endswith('.cram'), self.bam_or_cram_path
                index_path = self.bam_or_cram_path + '.crai'
                
        if isinstance(self.bam_or_cram_path, str):
            return b.read_input_group(
                base=self.bam_or_cram_path,
                index=index_path
            )
        else:
            return self.bam_or_cram_path
