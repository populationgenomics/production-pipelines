"""
wrapper concept for calling annotate_mito_coverage
"""


from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, copy_common_env
from cpg_workflows.batch import get_batch
from cpg_workflows.python_scripts import annotate_mito_coverage


# this logic probably be done in the relevant queue_jobs method
anno_job = get_batch().new_job(name="run mito annotate script")
# presumably use the same image as the rest of this run's hail jobs?
anno_job.image(get_config()['workflows']['driver_image'])

# might not be needed
copy_common_env(anno_job)

# authenticate google cloud credentials
authenticate_cloud_credentials_in_job(anno_job)

# command calls the script by file path
# consistent so long as _this_ job uses the same image (driver)
anno_job.command(
    f'python3 {annotate_mito_coverage.__file__} '
    '-i {input_tsv} '
    '-o {output_ht} '
    '-i {temp_dir} '   # could dump this for output_path(suffix, 'tmp')
)

# probably return the job to whatever called this

get_batch().run(wait=False)
