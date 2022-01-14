import os


PROJECT = 'fewgenomes'
BASE_BUCKET = 'gs://cpg-fewgenomes-test/unittest'

# Samples for joint calling
SAMPLES = [
    'CPG196519', 'CPG196527', 'CPG196535', 'CPG196543', 'CPG196550', 
    'CPG196568', 'CPG196576', 'CPG196584', 'CPG196592', 'CPG196600'
]
FULL_GVCF_BY_SID = {
    s: f'gs://cpg-fewgenomes-test/unittest/inputs/gvcf/{s}.g.vcf.gz'
    for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/gvcf/nagim/

SUBSET_GVCF_BY_SID = {
    s: f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/gvcf/{s}.g.vcf.gz'
    for s in SAMPLES
}


def setup_env(
    dataset: str = PROJECT,
    dataset_gcp_project: str = PROJECT,
    access_level: str = 'test',
):
    assert os.environ.get('HAIL_TOKEN')
    os.environ['DATASET'] = dataset
    os.environ['DATASET_GCP_PROJECT'] = dataset_gcp_project
    os.environ['ACCESS_LEVEL'] = access_level
    os.environ['HAIL_BUCKET'] = f'cpg-{dataset}-test-tmp/hail'
    os.environ['HAIL_BILLING_PROJECT'] = dataset
    os.environ['OUTPUT'] = 'unittests'


# Make sure setup_env() is called before the imports, because the env vars
# are required by the dataproc module on import time, and it's also
# imported by cpg_pipes
setup_env()
