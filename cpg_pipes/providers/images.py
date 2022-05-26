"""
Docker images used in the pipelines.
"""
import os


AVAILABLE = {
    'gatk': 'gatk:4.2.6.1',
    'bcftools': 'bcftools:1.10.2--h4f4756c_2',
    'sm-api': 'sm-api:4.0.0',
    'bwa': 'bwa:v0',
    'bwamem2': 'bwamem2:v0',
    'dragmap': 'dragmap:1.3.0',
    'samtools': 'picard_samtools:v0',
    'picard': 'picard_samtools:v0',
    'picard_samtools': 'picard_samtools:v0',
    'somalier': 'somalier:v0.2.15',
    'peddy': 'peddy:v0',
    'vep': 'vep:105',
    'verify-bam-id': 'verify-bam-id:1.0.1',
    'multiqc': 'multiqc:v1.12',
    'fastqc': 'fastqc:v0.11.9_cv8',
}


class Images:
    """
    Docker images used in the pipelines.
    """

    def __init__(self, prefix: str):
        self._prefix = prefix

    def get(self, name) -> str:
        """
        Get image URL by name
        """
        if name not in AVAILABLE:
            raise ValueError(
                f'Image {name} is not available. Known images: {AVAILABLE}'
            )
        return os.path.join(self._prefix, AVAILABLE[name])

    @staticmethod
    def driver_image() -> str:
        """
        Hail "driver" image.
        """
        return os.environ['CPG_DRIVER_IMAGE']
