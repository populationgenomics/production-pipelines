#!/usr/bin/env python3

import logging
from os.path import join
import pandas as pd
from analysis_runner import dataproc
from sample_metadata import (
    AnalysisApi,
    SampleApi,
)
from cpg_production_pipelines import utils
from cpg_production_pipelines.pipeline import Pipeline

sapi = SampleApi()
aapi = AnalysisApi()

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test-analysis/benchmark_combiner'

INPUT_PROJECT = 'fewgenomes'

pipe = Pipeline(
    analysis_project='fewgenomes',
    name='benchmark_combiner',
    output_version='v0',
    namespace='test',
    title='Benchmark GVCF combiner',
    smdb_check_existence=False,
)
pipe.find_samples(input_projects=[INPUT_PROJECT], namespace='test')
samples = pipe.get_all_samples()
gvcf_analysis_per_sid = pipe.db.find_analyses_by_sid(
    sample_ids=[s.id for s in samples],
    analysis_type='gvcf',
    project=INPUT_PROJECT,
)
df = pd.DataFrame(
    [
        {'s': s.id, 'gvcf': gvcf_analysis_per_sid[s.id]}
        for s in samples
        if s.id in gvcf_analysis_per_sid
    ]
)
logger.info(f'Found {len(df)} samples in {INPUT_PROJECT} with gvcfs')

for n_workers in [20, 40]:
    for n_samples in [100, 200, 400]:
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

pipe.run()
