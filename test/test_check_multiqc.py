"""
Test Hail Query functions.
"""

import toml
from pytest_mock import MockFixture
from cpg_utils import to_path

TOML = """
[workflow]
sequencing_type = 'genome'

[qc_thresholds.genome.min]
"MEDIAN_COVERAGE" = 10
"PCT_PF_READS_ALIGNED" = 80
[qc_thresholds.genome.max]
"FREEMIX" = 0.04
"PERCENT_DUPLICATION" = 25
"""


def test_check_multiqc(mocker: MockFixture, caplog):
    mocker.patch('cpg_utils.config.get_config', lambda: toml.loads(TOML))

    from cpg_workflows.python_scripts import check_multiqc

    data_dir = to_path(__file__).parent / 'data' / 'check_multiqc'
    check_multiqc.run(str(data_dir / 'validation_multiqc.json'))
    for expected_line in ['⭕ CPG243717|NA12878_KCCG: MEDIAN_COVERAGE=8.00<10.00']:
        matching_lines = [expected_line in msg for msg in caplog.messages]
        assert any(matching_lines)
