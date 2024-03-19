"""
Test Hail Query functions.
"""

from cpg_utils import to_path

from . import set_config

TOML = """
[workflow]
sequencing_type = 'genome'

[qc_thresholds.genome.min]
"MEDIAN_COVERAGE" = 10
"PCT_PF_READS_ALIGNED" = 0.80
[qc_thresholds.genome.max]
"FREEMIX" = 0.04
"PERCENT_DUPLICATION" = 25
"""


def test_check_multiqc(caplog, tmp_path):
    set_config(TOML, tmp_path / 'config.toml')

    from cpg_workflows.python_scripts import check_multiqc

    data_dir = to_path(__file__).parent / 'data' / 'check_multiqc'
    check_multiqc.run(str(data_dir / 'validation_multiqc.json'), send_to_slack=False)
    assert '❗ CPGeee|NA12878_KCCG: MEDIAN_COVERAGE=8.00<10.00' in caplog.messages
