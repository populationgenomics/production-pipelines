"""
Test Hail Query functions.
"""

import toml
from pytest_mock import MockFixture
from cpg_utils import to_path

TOML = """
[workflow]
sequencing_type = 'genome'
"""


def test_check_pedigree(mocker: MockFixture, caplog):
    mocker.patch('cpg_utils.config.get_config', lambda: toml.loads(TOML))

    from cpg_workflows.python_scripts import check_pedigree

    data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'
    check_pedigree.run(
        somalier_samples_fpath=str(data_dir / 'somalier-samples.tsv'),
        somalier_pairs_fpath=str(data_dir / 'somalier-pairs.tsv'),
        expected_ped_fpath=str(data_dir / 'samples.ped'),
    )
    for expected_line in [
        '4/51 PED samples with mismatching sex',
        'Sex inferred for 51/51 samples, matching for 47 samples.',
        'Found 4 sample pair(s) that are provided as unrelated, are inferred as related:',
        'Found 7 sample pair(s) that are provided as related, but inferred as unrelated:',
    ]:
        matching_lines = [expected_line in msg for msg in caplog.messages]
        assert any(matching_lines)
