"""
Reference files and indices used in bionformatics pipelines.
"""
from cpg_pipes.providers.refdata import RefData
from cpg_utils.hail_batch import reference_path


class CpgRefData(RefData):
    """
    Reference files in CPG.
    """
    
    def __init__(self):
        super().__init__(reference_path(suffix=''))
