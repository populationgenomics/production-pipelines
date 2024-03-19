#!/usr/bin/env python3

"""
Benchmarking Hail VDS combiner. Run with analysis runner:

analysis-runner --dataset thousand-genomes --access-level standard <script>
"""

import logging

import pandas as pd
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, output_path

from cpg_workflows import get_batch, get_cohort

# benchmark matrix:
N_WORKERS = [10, 30, 20, 40, 50]
N_SEQUENCING_GROUPS = [100, 200, 300, 400, 500]


def main():
    df = pd.DataFrame([{'s': s.id, 'gvcf': s.gvcf} for s in get_cohort().get_sequencing_groups() if s.gvcf.exists()])
    logging.info(
        f'Found {len(df)}/{len(get_cohort().get_sequencing_groups())} sequencing groups '
        f'in {get_config()["workflow"]["dataset"]} with GVCFs'
    )

    out_prefix = to_path(dataset_path('benchmark-combiner'))
    tmp_prefix = to_path(dataset_path('benchmark-combiner', category='tmp'))

    for n_workers in N_WORKERS:
        for n_sequencing_groups in N_SEQUENCING_GROUPS:
            label = f'nseqgroups-{n_sequencing_groups}-nworkers-{n_workers}'
            out_vds_path = out_prefix / label / 'combined.vds'
            sequencing_group_ids = get_cohort().get_sequencing_group_ids()[:n_sequencing_groups]

            from cpg_workflows.large_cohort.combiner import run
            from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

            dataproc_job(
                job_name=f'Combine {n_sequencing_groups} GVCFs on {n_workers} workers',
                function=run,
                function_path_args=dict(
                    out_vds_path=out_vds_path,
                    tmp_prefix=tmp_prefix,
                ),
                function_str_args=sequencing_group_ids,
                autoscaling_policy=(get_config()['hail'].get('dataproc', {}).get('combiner_autoscaling_policy')),
                num_workers=n_workers,
            )

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
