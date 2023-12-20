from cpg_utils.hail_batch import get_batch, output_path

from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
from cpg_workflows.large_cohort.sample_qc import run

dataproc_job(
    job_name='SampleQC',
    function=run,
    function_path_args=dict(
        vds_path='gs://cpg-bioheart-test/vds/3-1-3.vds',
        out_sample_qc_ht_path=output_path('sample-qc-out'),
        tmp_prefix=output_path('2023-12-18-Ick', 'tmp'),
    ),
)

get_batch().run(wait=False)
