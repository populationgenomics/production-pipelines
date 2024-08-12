"""
Generate new intervals from a MatrixTable based on data distribution
A wrapper for the generate_new_intervals.py script
"""

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch

# some hard coded paths
save_bed = 'gs://cpg-seqr-test/seqr_loader/new_intervals.bed'
input_mt = 'gs://cpg-seqr-test/seqr_loader/1d93d4819ca2d100c5108e4b0973915b5eb4a9_3905/AnnotateCohort/cohort.mt'

job = get_batch().new_bash_job('generate_intervals')
authenticate_cloud_credentials_in_job(job)
job.image(config_retrieve(['workflow', 'driver_image']))

# read this file in locally
meres = 'gs://cpg-common-test/references/hg38/v0/hg38.telomeresAndMergedCentromeres.interval_list'
meres_in = get_batch().read_input(meres)

job.command(f'new_intervals_from_mt --mt {input_mt} --out {job.output} --meres_file {meres_in} --intervals 1500')
get_batch().write_output(job.output, save_bed)
get_batch().run(wait=False)
