"""
Script to just populate analysis entries with status=completed in the 
sample-metadata DB. For back-populating old data; for new data, it should be
populated automatically with `workflow/status_provider="metamist"` set in config.
"""
from cloudpathlib.exceptions import OverwriteNewerCloudError
from cpg_utils import to_path

from cpg_pipes.inputs import get_cohort


cohort = get_cohort()

MOVE_DUPLICATE_METRICS = True
MOVE_CRAM_QC = False
DRY_RUN = False
REMOVE_METRICS = False

if REMOVE_METRICS:
    for dataset in cohort.get_datasets():
        for i, sample in enumerate(cohort.get_samples()):
            path = to_path(
                f'gs://cpg-{dataset.name}-main/exome/qc/verify_bamid/{sample.id}_verify_bamid.selfSM'
            )
            if path.exists():
                print(f'{i} removing {path}')
                path.unlink()
                assert not path.exists()

if MOVE_DUPLICATE_METRICS:
    for i, sample in enumerate(cohort.get_samples()):
        current_path = (
            sample.make_cram_path().path.parent
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
            print(f'{i} Warning: {current_path} does not exist')
        else:
            print(f'{i} cp {current_path}->{new_path}')
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
                print(f'{i} Warning: {old_path} does not exist')
            else:
                print(f'{i} Copy from {old_path} to {new_path}')
                if not DRY_RUN:
                    try:
                        old_path.copy(new_path)
                    except OverwriteNewerCloudError:
                        pass
