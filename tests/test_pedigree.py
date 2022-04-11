"""
Testing pedigree checks.
"""

import unittest

from cpg_pipes import to_path
from cpg_pipes.jobs.scripts import check_pedigree


class TestPedigree(unittest.TestCase):
    """
    Test pedigree checks.
    """

    def test_csv_provider(self):
        data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'

        check_pedigree.check_pedigree(
            somalier_samples_fpath=str(data_dir / 'somalier-samples.tsv'),
            somalier_pairs_fpath=str(data_dir / 'somalier-pairs.tsv'),
            sample_map_tsv_path=str(data_dir / 'sample-map.tsv'),
        )


if __name__ == '__main__':
    unittest.main()
