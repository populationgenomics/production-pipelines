"""
Test Hail Query functions.
"""
import shutil

import hail as hl
import pytest
import toml

from cpg_utils import to_path, Path
from cpg_utils.config import get_config, set_config_paths, update_dict
from cpg_utils.workflows.batch import get_batch
from cpg_utils.workflows.utils import timestamp
from cpg_utils.hail_batch import dataset_path, init_batch

from jobs.python_scripts import check_pedigree


def test_check_pedigree_good():
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
        matching_lines = [expected_line in msg for msg in result.messages]
        assert matching_lines


def test_check_pedigree_bad():
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
        matching_lines = [expected_line in msg for msg in result.messages]
        assert matching_lines
