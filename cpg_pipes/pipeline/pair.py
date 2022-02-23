"""
Pair of samples as a stage's Target
"""
from cpg_pipes.pipeline.target import Target
from cpg_pipes.pipeline.sample import Sample


class Pair(Target):
    """
    Pair of samples
    """

    @property
    def unique_id(self) -> str:
        return f'{self.s1.unique_id}:{self.s2.unique_id}'

    def __init__(self, s1: Sample, s2: Sample):
        super().__init__()
        self.s1 = s1
        self.s2 = s2
