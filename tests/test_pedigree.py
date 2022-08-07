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

    def test_check_pedigree_good(self):
        data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'
        result = check_pedigree.check_pedigree(
            somalier_samples_fpath=str(data_dir / 'output_samples_good.tsv'),
            somalier_pairs_fpath=str(data_dir / 'output_pairs_good.tsv'),
            expected_ped_fpath=str(data_dir / 'output_samples_good.tsv'),
        )
        for expected_line in [
            'Sex inferred for 28/28 samples, matching for all samples',
            'Inferred pedigree matches for all 28 samples',
        ]:
            self.assertTrue(any(expected_line in msg for msg in result.messages))

    def test_check_pedigree_bad(self):
        data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'
        result = check_pedigree.check_pedigree(
            somalier_samples_fpath=str(data_dir / 'output_samples_bad.tsv'),
            somalier_pairs_fpath=str(data_dir / 'output_pairs_bad.tsv'),
            expected_ped_fpath=str(data_dir / 'output_samples_bad.tsv'),
        )
        for expected_line in [
            '1/27 PED samples with mismatching sex',
            '2/27 samples with missing provided sex',
            '2/27 samples with failed inferred sex',
            'Sex inferred for 25/27 samples, matching for 24 samples',
            'Found 1 sample pair(s) that are provided as related, but inferred as unrelated:',
        ]:
            self.assertTrue(any(expected_line in msg for msg in result.messages))


if __name__ == '__main__':
    unittest.main()
