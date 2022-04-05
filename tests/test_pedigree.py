"""
Testing CSV inputs provider. CSV file must contain one required column "sample".
"""
import time
import unittest

from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.jobs.somalier import check_pedigree_job

try:
    from .utils import BASE_BUCKET, DATASET, SAMPLES
except ImportError:
    from utils import BASE_BUCKET, DATASET, SAMPLES  # type: ignore


class TestPedigree(unittest.TestCase):
    """
    Test Pedigree checks.
    """
    def setUp(self) -> None:
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.b = setup_batch(
            self.name, 
            billing_project=DATASET,
            hail_bucket=self.tmp_bucket
        )

    @property
    def tmp_bucket(self):
        return BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def test_check_pedigree(self):
        inputs_bucket = BASE_BUCKET / 'inputs' / 'check_pedigree/'
        check_pedigree_job(
            self.b,
            sample_map_file=self.b.read_input(str(inputs_bucket / 'sample-map.tsv')),
            samples_file=self.b.read_input(str(inputs_bucket / 'somalier-samples.tsv')),
            pairs_file=self.b.read_input(str(inputs_bucket / 'somalier-pairs.tsv')),
        )
        self.b.run()


if __name__ == '__main__':
    unittest.main()
