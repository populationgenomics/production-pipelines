"""
Common test utilities.
"""
import toml
from cpg_utils import Path, to_path
from cpg_utils.config import set_config_paths, update_dict
from cpg_utils.flows.filetypes import GvcfPath, FastqPair

DATASET = 'fewgenomes'
ACCESS_LEVEL = 'test'
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

SEQ_TYPE = 'genome'


def setup_env(tmp_bucket: Path, extraconf: dict | None = None):
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
        },
        'hail': {
            'billing_project': DATASET,
        },
        'elasticsearch': {
            'password_project_id': 'TEST',
            'password_secret_id': 'TEST',
        },
    }
    update_dict(conf, extraconf or {})
    config_path = tmp_bucket / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(conf, f)
    set_config_paths([str(config_path)])
