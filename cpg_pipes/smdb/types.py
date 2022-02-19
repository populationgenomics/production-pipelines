"""
Sample-metadata DB types. 

Re-defined in a separate module to decouple from the main smdb module, so 
decorators can use `@stage(analysis_type=AnalysisType.QC)`
"""
from collections import Set
from dataclasses import dataclass
from enum import Enum
from typing import Optional


class AnalysisType(Enum):
    """
    Corresponds to SMDB Analysis types:
https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L4-L11
    """
    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'


@dataclass
class Analysis:
    """
    Sample-metadata DB "Analysis" entry.

    See sample-metadata for more details: https://github.com/populationgenomics/sample-metadata
    """
    id: int
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]
    
    
class SmSequence:
    """
    Sample-metadata DB "Sequence" entry.

    See sample-metadata for more details: https://github.com/populationgenomics/sample-metadata
    """
    def __init__(self, id, sample_id, meta, smdb):
        self.id = id
        self.sample_id = sample_id
        self.meta = meta
        self.smdb = smdb
