from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.query_modules import seqr_loader

subset_j = get_batch().new_job('Subset cohort to dataset')
subset_j.image(config_retrieve(['workflow', 'driver_image']))
subset_j.cpu(4).memory('highmem')
sequencing_group_ids = config_retrieve(['workflow', 'sgids'])
mt_path = 'gs://cpg-seqr-test/seqr_loader/1d93d4819ca2d100c5108e4b0973915b5eb4a9_3905/AnnotateCohort/cohort.mt'
subset_mt_path = 'gs://cpg-seqr-test/subset_cohort.mt'

assert sequencing_group_ids
subset_j.command(
    query_command(
        seqr_loader,
        seqr_loader.subset_mt_to_samples.__name__,
        mt_path,
        sequencing_group_ids,
        subset_mt_path,
        setup_gcp=True,
    ),
)
get_batch().run(wait=False)
