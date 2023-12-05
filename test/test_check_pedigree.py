"""
Test Hail Query functions.
"""

from cpg_utils import to_path

from . import set_config

TOML = """
[workflow]
sequencing_type = 'genome'
"""


def test_check_pedigree(caplog, tmp_path):
    from cpg_workflows.python_scripts import check_pedigree

    set_config(TOML, tmp_path / 'config.toml')

    data_dir = to_path(__file__).parent / 'data' / 'check_pedigree'
    check_pedigree.run(
        somalier_samples_fpath=str(data_dir / 'somalier-samples.tsv'),
        somalier_pairs_fpath=str(data_dir / 'somalier-pairs.tsv'),
        expected_ped_fpath=str(data_dir / 'samples.ped'),
        send_to_slack=False,
    )
    for expected_line in [
        '4/51 PED samples with mismatching sex',
        'Sex inferred for 51/51 samples, matching for 47 samples.',
        'Found 1 sample pair(s) that are provided as unrelated, are inferred as related:',
        'Found 7 sample pair(s) that are provided as related, but inferred as unrelated:',
    ]:
        matching_lines = [expected_line in msg for msg in caplog.messages]
        assert any(matching_lines)
