"""
Script to just populate analysis entries with status=completed in the 
sample-metadata DB. For back-populating old data; for new data, it should be
populated automatically with `workflow/status_provider="smdb"` set in config.
"""

from cpg_utils.config import get_config
from cpg_utils.hail_batch import Namespace

from cpg_pipes.providers.cpg.inputs import CpgInputProvider
from cpg_pipes.providers.cpg.smdb import SMDB
from cpg_pipes.providers.cpg.status import CpgStatusReporter
from cpg_pipes.targets import Cohort
from cpg_pipes.utils import exists

access_level = get_config()['workflow']['access_level']
cohort = Cohort(
    analysis_dataset_name=get_config()['workflow']['dataset'],
    namespace=Namespace.from_access_level(access_level),
)
smdb = SMDB(cohort.analysis_dataset.name)
input_provider = CpgInputProvider(smdb)
input_provider.populate_cohort(
    cohort=cohort,
    dataset_names=get_config()['workflow'].get('datasets'),
    skip_samples=get_config()['workflow'].get('skip_samples'),
    only_samples=get_config()['workflow'].get('only_samples'),
    skip_datasets=get_config()['workflow'].get('skip_datasets'),
)

status = CpgStatusReporter(smdb)

for sample in cohort.get_samples():
    if (path := sample.get_cram_path().path).exist():
        status.create_analysis(
            str(path),
            analysis_type='cram',
            analysis_status='completed',
            target=sample,
            meta=sample.get_job_attrs() | dict(
                size=path.stat().st_size,
                sequencing_type=cohort.sequencing_type.value,
            ),
        )
    if (path := sample.get_gvcf_path().path).exist():
        status.create_analysis(
            str(path),
            analysis_type='gvcf',
            analysis_status='completed',
            target=sample,
            meta=sample.get_job_attrs() | dict(
                size=path.stat().st_size,
                sequencing_type=cohort.sequencing_type.value,
            ),
        )

h = cohort.alignment_inputs_hash()
path = cohort.analysis_dataset.prefix() / 'mt' / f'{h}.mt'
if exists(path):
    status.create_analysis(
        str(path),
        analysis_type='joint-calling',
        analysis_status='completed',
        target=cohort,
        meta=cohort.get_job_attrs() | dict(
            sequencing_type=cohort.sequencing_type.value,
        ),
    )
