"""
Converts paths into Hail Batch inputs: ResourceFile and ResourceGroup objects.
"""

from pathlib import Path
from typing import cast

import hailtop.batch as hb
from hailtop.batch import ResourceGroup

from cpg_pipes import ref_data
from cpg_pipes.pipeline.analysis import CramPath


def fasta_group(b: hb.Batch) -> ResourceGroup:
    """
    Returns fasta reference resource group
    """
    return b.read_input_group(**ref_data.REF_D)


class AlignmentInput:
    """
    Represents inputs for an alignment job, which can be a set of fastq files,
    or a CRAM or a BAM file with an index.
    """
    def __init__(
        self, 
        fqs1: list[str|hb.ResourceFile]|None = None,
        fqs2: list[str|hb.ResourceFile]|None = None,
        cram_path: CramPath|Path|str|hb.ResourceGroup|None = None,
        index_path: Path|str|None = None,
    ):
        self.fqs1 = fqs1
        self.fqs2 = fqs2
        self.cram_path: CramPath|hb.ResourceGroup|None = None
        if isinstance(cram_path, str|Path):
            self.cram_path = CramPath(cram_path, index_path=index_path)
        else:
            self.cram_path = cram_path

    def __repr__(self):
        return (
            f'AlignmentInput(' +
            (self.cram_path if self.is_bam_or_cram() else (self.fqs1, self.fqs2)) +
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
            elif any(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]), self
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

    def get_fqs1(self) -> list[str|hb.ResourceFile]:
        assert self.is_fastq()
        return cast(list, self.fqs1)

    def get_fqs2(self) -> list[str|hb.ResourceFile]:
        assert self.is_fastq()
        return cast(list, self.fqs2)

    def as_fq_inputs(self, b) -> tuple[list[hb.Resource], list[hb.Resource]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        assert self.is_fastq()
        self.fqs1 = cast(list, self.fqs1)
        self.fqs2 = cast(list, self.fqs2)
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

        if isinstance(self.cram_path, hb.ResourceGroup):
            return cast(hb.ResourceGroup, self.cram_path)

        self.cram_path = cast(CramPath, self.cram_path)
        return self.cram_path.resource_group(b)
