import os
import string
import time
from random import random, choices

from cpg_pipes import to_path
from cpg_pipes.types import GvcfPath, FastqPair

DATASET = 'fewgenomes'
BASE_BUCKET = to_path('gs://cpg-fewgenomes-test/unittest')

# Samples for joint calling
SAMPLES = [
    'CPG196519',
    'CPG196527',
    'CPG196535',
    'CPG196543',
    'CPG196550',
    'CPG196568',
    'CPG196576',
    'CPG196584',
    'CPG196592',
    'CPG196600',
]

FULL_CRAM_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/cram/{s}.cram') for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/cram
TOY_CRAM_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/toy/cram/{s}.cram') for s in SAMPLES
}
# CRAM above converted into FASTQ
TOY_FQ_BY_SID = {
    s: FastqPair(
        BASE_BUCKET / f'inputs/toy/fq/{s}_R1.fastq.gz',
        BASE_BUCKET / f'inputs/toy/fq/{s}_R2.fastq.gz',
    )
    for s in SAMPLES
}

FULL_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/gvcf
EXOME_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/exome/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}
CHR20_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/chr20/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}


def timestamp():
    """
    Generate a timestamp plus a short random string for guaranteed uniqueness.
    """
    rand_bit = ''.join(choices(string.ascii_uppercase + string.digits, k=3))
    return time.strftime('%Y-%m%d-%H%M') + rand_bit


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
    token = os.environ.get('HAIL_TOKEN')
    assert token
    os.environ['HAIL_TOKEN'] = token
    os.environ['CPG_DATASET'] = dataset
    os.environ['CPG_OUTPUT_PREFIX'] = 'unittests'
    os.environ['CPG_ACCESS_LEVEL'] = access_level
    os.environ['CPG_DATASET_GCP_PROJECT'] = dataset_gcp_project
    os.environ['HAIL_BILLING_PROJECT'] = dataset


setup_env()
