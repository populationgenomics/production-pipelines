#!/usr/bin/env python3
from cpg_workflows.metamist import get_metamist
from cpg_utils.config import set_config_paths
import os

TOML = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
sequencing_type = 'genome'
check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'
[storage.default]
analysis = "gs://cpg-fewgenomes-test-analysis"
default = "gs://cpg-fewgenomes-test"
tmp = "gs://cpg-fewgenomes-test-tmp"
upload = "gs://cpg-fewgenomes-test-upload"
web = "gs://cpg-fewgenomes-test-web"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"
[storage.vb-fewgenomes]
analysis = "gs://cpg-fewgenomes-test-analysis"
default = "gs://cpg-fewgenomes-test"
tmp = "gs://cpg-fewgenomes-test-tmp"
upload = "gs://cpg-fewgenomes-test-upload"
web = "gs://cpg-fewgenomes-test-web"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"
[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
"""

config_paths = [
    'configs/genome.toml',
]

workflow = 'large_cohort'
wfl_conf_path = f'configs/defaults/{workflow}.toml'
config_paths = os.environ['CPG_CONFIG_PATH'].split(',') + list(config_paths)
set_config_paths(config_paths[:1] + [str(wfl_conf_path)] + config_paths[1:])

sequencing_group_entries = get_metamist().get_sg_entries('fewgenomes')
print(sequencing_group_entries)
