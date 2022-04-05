#!/usr/bin/env python3

"""
Generate sample map to upload a dataset to seqr.
"""
import logging
import tempfile

import click
import pandas as pd
from typing import Optional

from cpg_pipes import Path, Namespace, to_path
from cpg_pipes.pipeline.cli_opts import choice_from_enum, val_to_enum
from cpg_pipes.providers import Cloud
from cpg_pipes.providers.cpg import SmdbInputProvider, CpgStorageProvider
from cpg_pipes.providers.cpg.smdb import SMDB
from cpg_pipes.targets import Dataset, Cohort

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.option('--dataset', 'datasets', multiple=True)
@click.option(
    '-n', '--namespace', 'namespace',
    type=choice_from_enum(Namespace),
    callback=val_to_enum(Namespace),
    help='The bucket namespace to write the results to',
)
@click.option(
    '--cloud',
    type=choice_from_enum(Cloud),
    callback=val_to_enum(Cloud),
    default=Cloud.GS.value,
    help='Cloud storage provider',
)
def main(datasets: list[str], namespace: Namespace, cloud: Cloud): 
    """
    Entry point.
    """
    analysis_dataset = 'seqr'
    input_provider = SmdbInputProvider(SMDB(analysis_dataset))
    cohort = input_provider.populate_cohort(
        cohort=Cohort(
            name=analysis_dataset,
            analysis_dataset_name=analysis_dataset,
            namespace=namespace,
            storage_provider=CpgStorageProvider(cloud)
        ),
        dataset_names=datasets,
    )

    tmp_dir = to_path(tempfile.mkdtemp())
    for dataset in cohort.get_datasets():
        _make_seqr_metadata_files(
            dataset=dataset,
            bucket=cohort.analysis_dataset.get_analysis_bucket(),
            local_dir=tmp_dir,
        )


def _make_seqr_metadata_files(
    dataset: Dataset, 
    bucket: Path,
    local_dir: Path, 
):
    """
    Create Seqr metadata files
    """
    samplemap_bucket_path = bucket / 'seqr' / f'{dataset.name}-sample-map.csv'
    igv_paths_path = local_dir / f'{dataset.name}-igv-paths.tsv'

    # Sample map
    df = pd.DataFrame({
        'cpg_id': s.id,
        'individual_id': s.external_id,
    } for s in dataset.get_samples())
    df.to_csv(str(samplemap_bucket_path), sep=',', index=False, header=False)

    # IGV
    df = pd.DataFrame({
        'individual_id': s.external_id,
        'cram_path': s.get_cram_path(),
        'cram_sample_id': s.id,
    } for s in dataset.get_samples() if s.get_cram_path())
    df.to_csv(str(igv_paths_path), sep='\t', index=False, header=False)

    logger.info(f'{dataset.name} sample map: {samplemap_bucket_path}')
    logger.info(f'{dataset.name} IGV paths: {igv_paths_path}')


if __name__ == '__main__':
    main()  # pylint: disable=E1120
