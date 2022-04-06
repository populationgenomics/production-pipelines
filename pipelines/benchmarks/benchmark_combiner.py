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

from cpg_pipes import utils, Namespace
from cpg_pipes.pipeline import create_pipeline

sapi = SampleApi()
aapi = AnalysisApi()

logger = logging.getLogger(__file__)

# INPUT_DATASET = 'fewgenomes'
INPUT_DATASET = 'tob-wgs'
NAMESPACE = Namespace.MAIN
BENCHMARK_BUCKET = f'gs://cpg-{INPUT_DATASET}-{NAMESPACE}-analysis/benchmark_combiner'


pipe = create_pipeline(
    analysis_dataset=INPUT_DATASET,
    name='benchmark_combiner',
    namespace=NAMESPACE,
    datasets=[INPUT_DATASET], 
)

df = pd.DataFrame(
    [
        {'s': s.id, 'gvcf': s.get_gvcf_path()} 
        for s in pipe.get_all_samples()
        if s.get_gvcf_path().exists()
    ]
)
logger.info(
    f'Found {len(df)}/{len(pipe.get_all_samples())} samples '
    f'in {INPUT_DATASET} with GVCFs'
)

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
