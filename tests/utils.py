"""
Common test utilities.
"""
import string
import time
from random import choices

import toml
from cpg_utils import Path
from cpg_utils.config import set_config_paths

from cpg_pipes import to_path, Namespace
from cpg_pipes.types import GvcfPath, FastqPair, SequencingType

DATASET = 'fewgenomes'
ACCESS_LEVEL = 'test'
Namespace = Namespace.TEST
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
    s: GvcfPath(BASE_BUCKET / f'inputs/full/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}  # from gs://cpg-thousand-genomes-test/gvcf
EXOME_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/exome/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}
EXOME_1PCT_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/exome1pct/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}
CHR20_GVCF_BY_SID = {
    s: GvcfPath(BASE_BUCKET / f'inputs/chr20/gvcf/{s}.g.vcf.gz') for s in SAMPLES
}

SEQ_TYPE = SequencingType.GENOME


def timestamp():
    """
    Generate a timestamp plus a short random string for guaranteed uniqueness.
    """
    rand_bit = ''.join(choices(string.ascii_uppercase + string.digits, k=3))
    return time.strftime('%Y-%m%d-%H%M') + rand_bit


def setup_env(timestamp_: str, tmp_bucket: Path):
    """Create config for tests"""
    conf = {
        'workflow': {
            'dataset': DATASET,
            'dataset_gcp_project': DATASET,
            'check_intermediates': False,
            'check_expected_outputs': False,
            'access_level': 'test',
            'realignment_shards_num': 4,
            'hc_intervals_num': 4,
            'jc_intervals_num': 4,
            'vep_intervals_num': 4,
            'version': timestamp_,
        },
        'hail': {
            'billing_project': DATASET,
            'bucket': str(tmp_bucket),
        },
        'elasticsearch': {
            'password': 'TEST',
        }
    }
    config_path = tmp_bucket / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(conf, f)
    set_config_paths([str(config_path)])
