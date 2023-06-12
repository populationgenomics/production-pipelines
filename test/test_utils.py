"""
utils tests
"""

from cpg_utils import to_path
from cpg_workflows.utils import missing_from_pre_collected


def test_exists_on_pre_collected():
    """
    checks that we find missing files between two sets
    """

    test = {to_path('a'), to_path('b')}
    known = {to_path('a')}

    assert missing_from_pre_collected(test, known) == to_path('b')


def test_exists_on_pre_collected_multi():
    """
    checks that we find one of the missing files from two sets
    """

    test = {to_path('a'), to_path('b'), to_path('c'), to_path('d')}
    known = {to_path('a')}

    assert missing_from_pre_collected(test, known) in {
        to_path('b'),
        to_path('c'),
        to_path('d')
    }


def test_exists_on_pre_collected_negative():
    """
    checks that we don't find a missing file
    """

    test = {to_path('a'), to_path('b')}
    known = {to_path('a'), to_path('b')}

    assert missing_from_pre_collected(test, known) is None
