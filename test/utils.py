import json
import os

from analysis_runner.constants import ANALYSIS_RUNNER_PROJECT_ID
from cpg_utils.cloud import read_secret

from cpg_pipes.pipeline import Sample

PROJECT = 'fewgenomes'
BASE_BUCKET = 'gs://cpg-fewgenomes-test/unittest'

# Samples for joint calling
SAMPLES = [
    Sample(id, id) for id in [
        'CPG196519', 'CPG196527', 'CPG196535', 'CPG196543', 'CPG196550', 
        'CPG196568', 'CPG196576', 'CPG196584', 'CPG196592', 'CPG196600'
    ]
]
FULL_GVCF_BY_SID = {
    s.id: f'gs://cpg-fewgenomes-test/unittest/inputs/gvcf/{s.id}.g.vcf.gz'
    for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/gvcf/nagim/

SUBSET_GVCF_BY_SID = {
    s.id: f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/gvcf/{s.id}.g.vcf.gz'
    for s in SAMPLES
}


def get_server_config() -> dict:
    """Get the server-config from the secret manager"""
    return json.loads(read_secret(ANALYSIS_RUNNER_PROJECT_ID, 'server-config'))


def setup_env(
    dataset: str = PROJECT,
    dataset_gcp_project: str = PROJECT,
    access_level: str = 'test',
):
    server_config = get_server_config()
    hail_token = server_config[dataset][f'{access_level}Token']
    assert hail_token, server_config

    os.environ['DATASET'] = dataset
    os.environ['DATASET_GCP_PROJECT'] = dataset_gcp_project
    os.environ['ACCESS_LEVEL'] = access_level
    os.environ['HAIL_BUCKET'] = f'cpg-{dataset}-test-tmp/hail'
    os.environ['HAIL_BILLING_PROJECT'] = dataset
    os.environ['HAIL_TOKEN'] = hail_token
    os.environ['OUTPUT'] = 'unittests'
