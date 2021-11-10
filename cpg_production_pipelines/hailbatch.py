from dataclasses import dataclass
from typing import Optional, Union, List, Tuple
import hailtop.batch as hb


@dataclass
class BamOrCramAlignmentInput:
    """
    Represents inputs for an alignment job of CRAM or BAM files
    """
    bam_or_cram_path: Optional[Union[str, hb.ResourceGroup]] = None
    index_path: Optional[str] = None

    def get_cram_input_group(self, b) -> hb.ResourceGroup:
        """
        Makes a ResourceGroup of bam/cram with accompanying index
        """
        assert self.bam_or_cram_path
        if isinstance(self.bam_or_cram_path, str):
            return b.read_input_group(
                base=self.bam_or_cram_path,
                index=self.index_path or self.bam_or_cram_path + '.crai'
            )
        else:
            return self.bam_or_cram_path


@dataclass
class FqAlignmentInput:
    """
    Represents inputs for an alignment job of fastq files
    """
    fqs1: Optional[List[Union[str, hb.ResourceFile]]] = None
    fqs2: Optional[List[Union[str, hb.ResourceFile]]] = None

    def get_fq_inputs(self, b) -> Tuple[List[hb.Resource], List[hb.Resource]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        if isinstance(self.fqs1[0], hb.Resource):
            files1 = self.fqs1
            files2 = self.fqs2
        else:
            files1 = [b.read_input(f1) for f1 in self.fqs1]
            files2 = [b.read_input(f1) for f1 in self.fqs2]
        return files1, files2


@dataclass
class AlignmentInput(BamOrCramAlignmentInput, FqAlignmentInput):
    pass
