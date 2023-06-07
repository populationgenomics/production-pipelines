"""
utils tests
"""

from cpg_utils import to_path
from cpg_workflows.utils import exists_on_pre_collected


def test_exists_on_pre_collected():
    """
    checks that we find missing files from two lists
    """

    test = {to_path('a'), to_path('b')}
    known = {to_path('a')}

    assert exists_on_pre_collected(test, known) == to_path('b')


def test_exists_on_pre_collected_negative():
    """
    checks that we don't find a missing file
    """

    test = {to_path('a'), to_path('b')}
    known = {to_path('a'), to_path('b')}

    assert exists_on_pre_collected(test, known) is None
