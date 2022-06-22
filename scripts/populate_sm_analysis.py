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

POPULATE_SAMPLES = True
POPULATE_JOINT_CALL = False
POPULATE_ES_INDEX = False


if POPULATE_SAMPLES:
    for sample in cohort.get_samples():
        if (path := sample.get_cram_path().path).exists():
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
        if (path := sample.get_gvcf_path().path).exists():
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

if POPULATE_JOINT_CALL:
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

if POPULATE_ES_INDEX:
    for name in [
        'acute-care-genome-2022_0620_1843_l4h8u',
        'ravenscroft-arch-genome-2022_0618_1137_4qfyn',
        'circa-genome-2022_0618_1137_4qfyn',
        'ohmr3-mendelian-genome-2022_0618_1137_4qfyn',
        # 'validation-genome-2022_0618_1137_4qfyn',
        'mito-disease-genome-2022_0618_1137_4qfyn',
        'perth-neuro-genome-2022_0618_1137_4qfyn',
        'ohmr4-epilepsy-genome-2022_0618_1137_4qfyn',
        'hereditary-neuro-genome-2022_0618_1137_4qfyn',
        'ravenscroft-rdstudy-genome-2022_0618_1137_4qfyn',
        'heartkids-genome-2022_0618_1137_4qfyn',
    ]:
        ds_name = name.split('-genome-')[0]
        print(f'Adding {ds_name}')
        dataset = cohort.create_dataset(ds_name)
        status.create_analysis(
            str(name),
            analysis_type='es-index',
            analysis_status='completed',
            target=dataset,
            meta=dataset.get_job_attrs() | dict(
                sequencing_type=cohort.sequencing_type.value,
            ),
        )
