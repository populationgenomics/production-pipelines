#!/usr/bin/env python3

"""
Benchmarking VCF combiner.
"""

import logging
import sys
from os.path import join
import pandas as pd
from analysis_runner import dataproc
from sample_metadata.apis import (
    AnalysisApi,
    SampleApi,
)
from sample_metadata.models import (
    SampleUpdateModel,
)

from cpg_pipes import utils
from cpg_pipes.storage import Namespace
from cpg_pipes.pipeline.analysis import AnalysisType
from cpg_pipes.pipeline.pipeline import Pipeline

sapi = SampleApi()
aapi = AnalysisApi()

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

# INPUT_DATASET = 'fewgenomes'
INPUT_DATASET = 'tob-wgs'
NAMESPACE = Namespace.MAIN
BENCHMARK_BUCKET = f'gs://cpg-{INPUT_DATASET}-{NAMESPACE}-analysis/benchmark_combiner'


pipe = Pipeline(
    analysis_dataset=INPUT_DATASET,
    name='benchmark_combiner',
    output_version='v0',
    namespace=NAMESPACE,
    description='Benchmark GVCF combiner',
    check_smdb_seq=False,
    input_datasets=[INPUT_DATASET], 
)

samples = pipe.get_all_samples()
sample_ids = [s.id for s in samples]

gvcf_analysis_per_sid = pipe.db.find_analyses_by_sid(
    sample_ids=sample_ids,
    analysis_type=AnalysisType('gvcf'),
    dataset=INPUT_DATASET,
)
gvcf_by_sid = dict()
if gvcf_analysis_per_sid:
    logger.info(
        f'Found {len(gvcf_analysis_per_sid)} GVCF analysis enteies '
        f'in {INPUT_DATASET}'
    )
    gvcf_by_sid = {sid: a.output for sid, a in gvcf_analysis_per_sid.items()}
if not gvcf_analysis_per_sid:
    logger.info(
        f'Found no GVCF analysis for {INPUT_DATASET}, '
        f'looking for GVCFs in sequence meta'
    )
    seqs = pipe.db.seqapi.get_sequences_by_sample_ids(request_body=sample_ids)
    gvcf_seq_per_sid = {
        seq['sample_id']: seq for seq in seqs if seq['meta'].get('gvcf')
    }
    if not gvcf_seq_per_sid:
        logger.error('Could not find GVCF in sequence meta either')
        sys.exit(1)
    logger.info(
        f'Found {len(gvcf_seq_per_sid)} GVCF sequence entries '
        f'in {INPUT_DATASET}'
    )
    gvcf_by_sid = {
        sid: seq['meta']['gvcf'][0]['location'] 
        for sid, seq in gvcf_seq_per_sid.items()
    }
    for gvcf_path, sid in gvcf_by_sid.items():
        if gvcf_path.endswith('.haplotypeCalls.er.raw.g.vcf.gz'):
            pipe.db.sapi.update_sample(sid, SampleUpdateModel)
        
if not gvcf_by_sid:
    logger.error('No GVCFs found, exiting')
    sys.exit(1)
logger.info(f'Found {len(gvcf_by_sid)} GVCFs')

df = pd.DataFrame(
    [
        {'s': s.id, 'gvcf': gvcf_by_sid[s.id]}
        for s in samples
        if s.id in gvcf_by_sid
    ]
)
logger.info(f'Found {len(df)} samples in {INPUT_DATASET} with gvcfs')

for n_workers in [10, 30, 20, 40, 50]:
    for n_samples in [100, 200, 300, 400, 500]:
        label = f'nsamples-{n_samples}-nworkers-{n_workers}'
        out_mt_path = join(BENCHMARK_BUCKET, label, 'combined.mt')

        meta_csv_path = join(BENCHMARK_BUCKET, label, 'meta.csv')
        subset_df = df.sample(n=n_samples)
        logger.info(f'Subset dataframe to {len(subset_df)} samples')
        subset_df.to_csv(meta_csv_path, index=False, sep='\t', na_rep='NA')

        combiner_job = dataproc.hail_dataproc_job(
            pipe.b,
            f'{utils.QUERY_SCRIPTS_DIR}/combine_gvcfs.py '
            f'--meta-csv {meta_csv_path} '
            f'--out-mt {out_mt_path} '
            f'--bucket {BENCHMARK_BUCKET}/work',
            max_age='8h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=n_workers,
            job_name=f'Combine {n_samples} GVCFs on {n_workers} workers',
        )

pipe.submit_batch()
