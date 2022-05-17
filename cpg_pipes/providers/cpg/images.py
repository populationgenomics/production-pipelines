"""
Docker images used in the pipelines.
"""
from cpg_utils.hail_batch import image_path

from cpg_pipes.providers.images import Images


class CpgImages(Images):
    """
    Docker images used in the pipelines.
    """
    def __init__(self):
        super().__init__(image_path(''))
