#!/usr/bin/env python3

"""
Benchmarking Hail VDS combiner. Run with analysis runner:

analysis-runner --dataset thousand-genomes --access-level standard <script>
"""

import logging

import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, dataset_path
from cpg_workflows import get_cohort, get_batch


# benchmark matrix:
N_WORKERS = [10, 30, 20, 40, 50]
N_SAMPLES = [100, 200, 300, 400, 500]


def main():
    df = pd.DataFrame(
        [
            {'s': s.id, 'gvcf': s.make_gvcf_path()}
            for s in get_cohort().get_samples()
            if s.make_gvcf_path().exists()
        ]
    )
    logging.info(
        f'Found {len(df)}/{len(get_cohort().get_samples())} samples '
        f'in {get_config()["workflow"]["dataset"]} with GVCFs'
    )

    out_prefix = to_path(dataset_path('benchmark-combiner'))
    tmp_prefix = to_path(dataset_path('benchmark-combiner', category='tmp'))

    for n_workers in N_WORKERS:
        for n_samples in N_SAMPLES:
            label = f'nsamples-{n_samples}-nworkers-{n_workers}'
            out_vds_path = out_prefix / label / 'combined.vds'
            sample_ids = get_cohort().get_sample_ids()[:n_samples]

            from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
            from cpg_workflows.large_cohort.combiner import run

            dataproc_job(
                job_name=f'Combine {n_samples} GVCFs on {n_workers} workers',
                function=run,
                function_path_args=dict(
                    out_vds_path=out_vds_path,
                    tmp_prefix=tmp_prefix,
                ),
                function_str_args=sample_ids,
                autoscaling_policy=(
                    get_config()['hail']
                    .get('dataproc', {})
                    .get('combiner_autoscaling_policy')
                ),
                num_workers=n_workers,
            )

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
