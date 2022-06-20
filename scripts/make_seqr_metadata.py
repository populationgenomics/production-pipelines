#!/usr/bin/env python3

"""
Generate sample map to upload a dataset to Seqr
"""
import logging
import tempfile

import click
import pandas as pd
from cpg_utils.config import get_config

from cpg_pipes import Path, Namespace, to_path
from cpg_pipes.providers.cpg.inputs import CpgInputProvider
from cpg_pipes.providers.cpg.smdb import SMDB
from cpg_pipes.targets import Dataset, Cohort

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--use-participant-id/--use-external-id',
    'use_external_id',
    default=False,
    is_flag=True,
)
def main(
    use_external_id: bool = False,
):
    """
    Generate sample map to upload a dataset to Seqr
    """
    input_provider = CpgInputProvider(SMDB())
    cohort = Cohort(
        analysis_dataset_name='seqr',
        namespace=Namespace.from_access_level(get_config()['workflow']['access_level']),
    )
    input_provider.populate_cohort(
        cohort=cohort,
        dataset_names=get_config()['workflow']['datasets'],
    )

    tmp_dir = to_path(tempfile.mkdtemp())
    for dataset in cohort.get_datasets():
        _make_seqr_metadata_files(
            dataset=dataset,
            bucket=cohort.analysis_dataset.prefix(),
            local_dir=tmp_dir,
            use_external_id=use_external_id,
        )


def _make_seqr_metadata_files(
    dataset: Dataset,
    bucket: Path,
    local_dir: Path,
    use_external_id: bool,
):
    """
    Create Seqr metadata files
    """
    samplemap_bucket_path = bucket / 'sample-maps' / f'{dataset.name}-sample-map.csv'
    igv_paths_path = local_dir / f'{dataset.name}-igv-paths.tsv'

    # Sample map
    df = pd.DataFrame(
        {
            'cpg_id': s.id,
            'individual_id': s.external_id if use_external_id else s.participant_id,
        }
        for s in dataset.get_samples()
    )
    with samplemap_bucket_path.open('w') as fh:
        df.to_csv(fh, sep=',', index=False, header=False)

    # IGV
    df = pd.DataFrame(
        {
            'individual_id': s.external_id if use_external_id else s.participant_id,
            'cram_path': s.get_cram_path(),
            'cram_sample_id': s.id,
        }
        for s in dataset.get_samples()
        if s.get_cram_path()
    )
    with igv_paths_path.open('w') as fh:
        df.to_csv(fh, sep='\t', index=False, header=False)

    logger.info(f'{dataset.name} sample map: {samplemap_bucket_path}')
    logger.info(f'{dataset.name} IGV paths: {igv_paths_path}')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
