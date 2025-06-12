"""
A couple of quick validation tests
"""

import pytest

from cpg_workflows.stages.happy_validation import find_hash_from_path


def test_hash_finder():
    # terminal component of the path
    test_input = 'gs://a-bucket/directory/l0r3MIpsum_1234'
    assert find_hash_from_path(test_input) == 'l0r3MIpsum_1234'

    # followed by a file name
    test_input_2 = 'gs://a-bucket/directory/l0r3MIpsum_1234/file.ext'
    assert find_hash_from_path(test_input_2) == 'l0r3MIpsum_1234'


def test_hash_finder_fails():
    """"""
    with pytest.raises(ValueError):
        find_hash_from_path('')
