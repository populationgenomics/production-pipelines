"""
Script to just populate analysis entries with status=completed in the 
sample-metadata DB. For back-populating old data; for new data, it should be
populated automatically with `workflow/status_provider="smdb"` set in config.
"""
from cloudpathlib.exceptions import OverwriteNewerCloudError
from cpg_utils.config import get_config
from cpg_utils.hail_batch import Namespace

from cpg_pipes.providers.cpg.inputs import CpgInputProvider
from cpg_pipes.providers.cpg.smdb import SMDB
from cpg_pipes.targets import Cohort
from cpg_pipes.types import SequencingType

access_level = get_config()['workflow']['access_level']
cohort = Cohort(
    analysis_dataset_name=get_config()['workflow']['dataset'],
    namespace=Namespace.from_access_level(access_level),
    sequencing_type=SequencingType.parse(get_config()['workflow']['sequencing_type']),
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

MOVE_DUPLICATE_METRICS = False
MOVE_CRAM_QC = True
DRY_RUN = False


if MOVE_DUPLICATE_METRICS:
    for i, sample in enumerate(cohort.get_samples()):
        current_path = (
            sample.get_cram_path().path.parent
            / 'duplicate-metrics'
            / f'{sample.id}-duplicate-metrics.csv'
        )
        new_path = (
            sample.dataset.prefix()
            / 'qc'
            / 'markduplicates_metrics'
            / (sample.id + '.markduplicates-metrics')
        )
        if not current_path.exists():
            print(f'Warning: {current_path} does not exist')
        else:
            current_path.copy(new_path)

if MOVE_CRAM_QC:
    for i, sample in enumerate(cohort.get_samples()):
        old_path_d = {
            key: sample.dataset.prefix() / 'qc' / f'{sample.id}{suf}'
            for key, suf in {
                'samtools_stats': '_samtools_stats.txt',
                'picard_wgs_metrics': '_picard_wgs_metrics.csv',
                'verify_bamid': '_verify_bamid.selfSM',
            }.items()
        }
        new_path_d = {
            key: sample.dataset.prefix() / 'qc' / key / f'{sample.id}{suf}'
            for key, suf in {
                'samtools_stats': '_samtools_stats.txt',
                'picard_wgs_metrics': '_picard_wgs_metrics.csv',
                'verify_bamid': '_verify_bamid.selfSM',
            }.items()
        }
        for key in old_path_d.keys():
            old_path = old_path_d[key]
            new_path = new_path_d[key]
            if not old_path.exists():
                print(f'Warning: {old_path} does not exist')
            else:
                print(f'Copy from {old_path} to {new_path}')
                if not DRY_RUN:
                    try:
                        old_path.copy(new_path)
                    except OverwriteNewerCloudError:
                        pass
