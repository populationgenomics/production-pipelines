import os
import string
from random import choices
import time
from cpg_utils import to_path


def results_prefix():
    """
    Output directory for the test results.
    """
    path = (
        to_path(__file__).parent
        / 'results'
        / os.getenv(
            'TEST_TIMESTAMP',
            # Generate a timestamp string. Don't import `timestamp` from `cpg_utils`,
            # as it would break mocking in some tests.
            time.strftime('%Y_%m%d_%H%M')
            + '_'
            + ''.join(choices(string.ascii_uppercase + string.digits, k=5)),
        )
    ).absolute()
    path.mkdir(parents=True, exist_ok=True)
    return path


def update_dict(d1: dict, d2: dict) -> None:
    """
    Helper functon to merge one dict into another. We don't import `update_dict` from
    `cpg_utils`, as it would break mocking in some tests.
    """
    for k, v2 in d2.items():
        v1 = d1.get(k)
        if isinstance(v1, dict) and isinstance(v2, dict):
            update_dict(v1, v2)
        else:
            d1[k] = v2
