"""
Align RNA-seq reads to the genome using STAR.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.hail_batch import command
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    FastqPair,
    FastqPairs,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
import re


class STAR:
    """
    Construct a STAR command for aligning FASTQs.
    """

    def __init__(
            self,
            input_fastq_pair: FastqPair,
            sample_name: str,
            genome_dir: str | Path,
            nthreads: int,
            read_group: dict[str, str] | None = None,
            bamout: bool = True,
            sort: bool = True,
            stdout: bool = True,
    ):
        self.command = ['STAR']
        # Create read group line
        self.read_group = read_group or {}
        if 'ID' not in self.read_group:
            self.read_group['ID'] = sample_name
        if 'SM' not in self.read_group:
            self.read_group['SM'] = sample_name
        self.read_group_line = self._create_read_group_line(self.read_group)
        # Check genome_dir exists
        if not Path(genome_dir).exists():
            raise FileNotFoundError(f'Genome directory {genome_dir} does not exist')
        # Check FASTQs exist
        for fastq in input_fastq_pair:
            if not Path(fastq).exists():
                raise FileNotFoundError(f'FASTQ {fastq} does not exist')
        # Create outSAMtype and outStd strings
        if bamout and sort:
            self.outSAMtype = 'BAM SortedByCoordinate'
            self.outStd = 'BAM_SortedByCoordinate'
        elif bamout:
            self.outSAMtype = 'BAM Unsorted'
            self.outStd = 'BAM_Unsorted'
        elif sort:
            self.outSAMtype = 'SAM SortedByCoordinate'
            self.outStd = 'SAM'
        else:
            self.outSAMtype = 'SAM Unsorted'
            self.outStd = 'SAM'
        if not stdout:
            self.outStd = 'Log'
        # Create command
        self.command.extend([
            '--runThreadN', str(nthreads),
            '--genomeDir', str(genome_dir),
            '--outSAMtype', self.outSAMtype,
            '--outStd', self.outStd,
            '--outSAMattrRGline', self.read_group_line,
            '--readFilesIn', input_fastq_pair.r1, input_fastq_pair.r2,
        ])
        

    def __str__(self):
        return ' '.join(self.command)
    
    def __repr__(self):
        return self.__str__()

    def _create_read_group_line(self) -> str:
        """
        Create a read group line for the STAR command.
        """
        read_group_line = '@RG'
        for k, v in self.read_group.items():
            v_quoted = re.sub(r'(^.*\s.*$)', r'"\1"', v)
            read_group_line += f' {k}:{v_quoted}'
        return read_group_line
    

def align():
    pass

def align_fq_pair():
    pass

def merge_bams():
    pass

def sort_bam():
    pass