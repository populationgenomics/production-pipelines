import os

from cloudpathlib import CloudPath

from cpg_pipes.pipeline.analysis import GvcfPath

DATASET = 'fewgenomes'
BASE_BUCKET = CloudPath('gs://cpg-fewgenomes-test/unittest')

# Samples for joint calling
SAMPLES = [
    'CPG196519', 'CPG196527', 'CPG196535', 'CPG196543', 'CPG196550', 
    'CPG196568', 'CPG196576', 'CPG196584', 'CPG196592', 'CPG196600'
]
FULL_GVCF_BY_SID = {
    s: GvcfPath(f'gs://cpg-fewgenomes-test/unittest/inputs/gvcf/{s}.g.vcf.gz')
    for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/gvcf/nagim/

SUBSET_GVCF_BY_SID = {
    s: GvcfPath(f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/gvcf/{s}.g.vcf.gz')
    for s in SAMPLES
}


def setup_env(
    dataset: str = DATASET,
    dataset_gcp_project: str = DATASET,
    access_level: str = 'test',
):
    """
    Make sure that setup_env() is called before importing cpg_pipes, because 
    cpg_pipes imports analysis_runner.dataproc, which in turns requires
    DATASET_GCP_PROJECT to be set on import time. Also, unittest doesn't pick
    exported environment variables with EnvFile, so have to do it here.
    """
    assert os.environ.get('HAIL_TOKEN')
    os.environ['CPG_DATASET'] = dataset
    os.environ['CPG_DATASET_GCP_PROJECT'] = dataset_gcp_project
    os.environ['CPG_ACCESS_LEVEL'] = access_level
    os.environ['HAIL_BUCKET'] = f'cpg-{dataset}-test-tmp/hail'
    os.environ['HAIL_BILLING_PROJECT'] = dataset
    os.environ['CPG_OUTPUT_PREFIX'] = 'unittests'


setup_env()
