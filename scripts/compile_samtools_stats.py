#!/usr/bin/env python3

"""
A script to read the samtools stats files and extract a given statistic for any number of samples
"""

import csv
import logging
import sys
from io import StringIO

import click
import tqdm
from google.cloud import storage

from cpg_utils import Path, to_path
from metamist.apis import ProjectApi
from metamist.graphql import gql, query

handler = logging.StreamHandler()
formatter = logging.Formatter(
    fmt='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
handler.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(handler)
logger.setLevel(logging.INFO)
logger.propagate = False

proj_api = ProjectApi()


def get_datasets(seqr_datasets_only: bool) -> list[str]:
    """
    Get the datasets from the Metamist Project API.
    Use the provided argument to filter datasets, or return all datasets if no argument is provided.
    """
    if seqr_datasets_only:
        # If wanting only RD datasets, fetch using the get_seqr_projects endpoint
        return [
            dataset['name']
            for dataset in proj_api.get_seqr_projects()
            if 'test' not in dataset['name'] and 'training' not in dataset['name']
        ]

    # If no argument is provided, get all datasets
    return [
        dataset['name']
        for dataset in proj_api.get_my_projects()
        if 'test' not in dataset['name'] and 'training' not in dataset['name']
    ]


SGS_BY_DATASET_QUERY = gql(
    """
    query MyQuery( $projects: [String!] ) {
      sequencingGroups(
        project: {in_: $projects}
        type: {eq: "genome"}
        technology: {eq: "short-read"}
        ) {
            id
            sample {
                project {
                    name
                }
            }
        }
    }
    """,
)


def get_sequencing_groups_by_dataset(datasets: list[str], excluded_sgs: set[str]) -> dict[str, list[str]]:
    """
    Get the sequencing groups for each dataset
    """
    sgs_query_results = query(
        SGS_BY_DATASET_QUERY,
        variables={"projects": datasets},
    )["sequencingGroups"]

    sgs_by_dataset: dict[str, list[str]] = {}
    for sg in sgs_query_results:
        if sg['id'] in excluded_sgs:
            continue
        dataset_name = sg['sample']['project']['name']
        if dataset_name not in sgs_by_dataset:
            sgs_by_dataset[dataset_name] = []
        sgs_by_dataset[dataset_name].append(sg['id'])

    return sgs_by_dataset


def make_samtools_stats_filepath(dataset: str, sequencing_group_id: str) -> Path:
    """
    Make the samtools stats file path for a given dataset and sequencing group ID
    """
    return to_path(
        f"gs://cpg-{dataset}-main/qc/samtools_stats/{sequencing_group_id}.samtools-stats",
    )


def get_samtools_stats_filepaths(
    dataset_sgs: dict[str, list[str]],
) -> dict[str, Path]:
    """
    Get the samtools stats filepath for each sequencing group in each dataset
    """
    samtools_stats_filepaths_by_sg: dict[str, Path] = {}
    for dataset, sgs in dataset_sgs.items():
        for sg_id in sgs:
            samtools_stats_filepaths_by_sg[sg_id] = make_samtools_stats_filepath(dataset, sg_id)

    return samtools_stats_filepaths_by_sg


def get_samtools_stats_from_file(filepath: Path) -> dict[str, str]:
    """
    Read a samtools stats file and return the summary statistics
    """
    stats = {}
    with filepath.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith('FFQ'):
                # Stop reading, as the summary statistics are above this line
                break

            if not line.startswith('SN\t'):
                continue

            parts = line.split('\t')
            if len(parts) >= 3:
                field_name = parts[1].rstrip(':').strip()
                value = parts[2].strip()

                # Remove any trailing comments
                if '#' in value:
                    value = value.split('#')[0].strip()

                stats[field_name] = value

    return stats


def get_samtools_stats_for_sgs(
    dataset_sgs: dict[str, list[str]],
    samtools_stats_filepaths_by_sg: dict[str, Path],
) -> dict[str, dict[str, dict[str, str]]]:
    """
    Read the samtools stats files for each sequencing group and return the statistics

    Returns a dictionary with dataset names as keys, and another dictionary as values,
    where the inner dictionary has sequencing group IDs as keys and their stats as values.
    {'dataset_name': {'sg_id': {'stat_name': 'value', ...}, ...}, ...}
    """
    datasets_sg_stats: dict[str, dict[str, dict[str, str]]] = {}
    for dataset, sgs in dataset_sgs.items():
        datasets_sg_stats[dataset] = {}
        for sg_id in tqdm.tqdm(sgs, desc=f'{dataset} :: Reading samtools stats files...'):
            filepath = samtools_stats_filepaths_by_sg[sg_id]
            if not filepath.exists():
                logger.warning(f'Samtools stats file not found: {filepath}')
                continue

            stats = get_samtools_stats_from_file(filepath)
            datasets_sg_stats[dataset][sg_id] = stats

    return datasets_sg_stats


def get_fieldnames_from_stats(sg_stats: dict[str, str]) -> list[str]:
    """
    Get the fieldnames from the stats dictionary.
    This will return a list of all unique fieldnames across all stats.
    """
    fieldnames = ['dataset', 'sg_id']
    for stat in sg_stats.keys():
        if stat not in fieldnames:
            fieldnames.append(stat)
    return fieldnames


def get_existing_data(
    output_file: str,
    bucket_name: str | None,
    fieldnames: list[str],
) -> list[str]:
    """
    Read an existing output file and return its contents as a list of lines.
    """
    exists = to_path(output_file).exists()
    if exists:
        logger.info(f'Output file {output_file} already exists. Appending to it.')
        # Download the existing file and read its contents, then append the new stats to it and write it back
        if output_file.startswith('gs://'):
            # If the output file is in GCS, download it to a temporary local file
            storage_client = storage.Client()
            bucket = storage_client.bucket(bucket_name)
            blob = bucket.blob(output_file.removeprefix(f'gs://{bucket_name}/'))
            temp_file = to_path('/tmp/temp_samtools_stats.tsv')
            blob.download_to_filename(temp_file)
            with open(temp_file, 'r', encoding='utf-8') as f:
                logger.info(f'Reading existing data from {temp_file}')
                reader = csv.DictReader(f, delimiter='\t')
                if reader.fieldnames != fieldnames:
                    logger.warning(
                        f'Fieldnames in existing file {temp_file} do not match expected fieldnames. Expected: {fieldnames}, Found: {reader.fieldnames}',
                    )
                existing_data = f.readlines()
            temp_file.unlink()  # Clean up the temporary file
        else:
            # If the output file is local, read it directly
            logger.info(f'Reading existing data from {output_file}')
            # Open the file and read its contents, skipping the header line
            with open(output_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                if reader.fieldnames != fieldnames:
                    logger.warning(
                        f'Fieldnames in existing file {output_file} do not match expected fieldnames. Expected: {fieldnames}, Found: {reader.fieldnames}',
                    )
                # Read the existing data to avoid duplicating headers
                existing_data = f.readlines()
        logger.info(f'Found {len(existing_data)} existing lines in {output_file}')
        return existing_data
    else:
        logger.info(f'Output file {output_file} does not exist. Creating it.')
        return []


def write_all_samtools_stats_to_file(
    datasets_sg_stats: dict[str, dict[str, dict[str, str]]],
    bucket_name: str | None,
    output_file: str,
):
    """
    Write the samtools stats to a file in TSV format, or append to an existing file.
    """
    # Get the fieldnames from the first dataset's first sequencing group's stats
    # All datasets should have the same fieldnames, so we can just take the first one
    stats_dict = next(iter(datasets_sg_stats.values())).get(
        next(iter(next(iter(datasets_sg_stats.values())).keys())),
        {},
    )
    if not stats_dict:
        logger.error('No statistics found in the provided datasets. Exiting.')
        sys.exit(1)
    fieldnames = get_fieldnames_from_stats(stats_dict)

    existing_data = get_existing_data(
        output_file,
        bucket_name,
        fieldnames,
    )

    # Write the data to a buffer
    buffer = StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=fieldnames,
        delimiter='\t',
    )
    writer.writeheader()
    # Write existing data to the buffer
    for line in existing_data:
        writer.writerow(dict(zip(fieldnames, line.strip().split('\t'))))
    # Write the new stats to the buffer
    for dataset, sg_stats in datasets_sg_stats.items():
        for sg_id, stats in sg_stats.items():
            writer.writerow(
                {
                    'dataset': dataset,
                    'sg_id': sg_id,
                    **stats,
                },
            )

    # Save to a local file or GCS blob based on the output_file path
    if not bucket_name:
        # If no bucket name is provided, save to local file
        with open(output_file, 'w') as f:
            f.write(buffer.getvalue())
            buffer.close()
            logger.info(f'Stats table written to {output_file}')
        return

    # If a bucket name is provided, save to GCS
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(output_file.removeprefix(f'gs://{bucket_name}/'))
    blob.upload_from_string(buffer.getvalue(), content_type='text/tab-separated-values')
    buffer.close()
    logger.info(f'Stats table written to gs://{bucket_name}/{blob.name}')
    return


def read_existing_outputs_sgs(
    output_file: Path,
) -> set[str]:
    """
    Read an existing output file and return the SGs within
    """
    dataset_sgs: dict[str, set[str]] = {}
    with output_file.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['dataset'] not in dataset_sgs:
                dataset_sgs[row['dataset']] = set()
            dataset_sgs[row['dataset']].add(row['sg_id'])

    logger.info(
        f'Found {len(dataset_sgs)} datasets with {sum(len(sgs) for sgs in dataset_sgs.values())} total sequencing groups in existing output file: {output_file}',
    )

    all_sgs = set()
    for dataset, sgs in dataset_sgs.items():
        logger.info(f'    {dataset} :: {len(sgs)} sequencing groups')
        all_sgs.update(sgs)

    return all_sgs


@click.command()
@click.option('-o', '--output-file', type=str)
@click.option('--rd', is_flag=True, help='If set (and --datasets is not set) use RD datasets. Else use all datasets.')
@click.option(
    '-d',
    '--datasets',
    type=str,
    multiple=True,
    help='Filter to specific datasets. If not provided, all datasets will be used.',
)
def main(output_file: str, rd: bool, datasets: list[str]) -> None:
    """
    Main function to run the script
    """
    if not output_file:
        logger.error('Output file must be specified with -o or --output-file')
        sys.exit(1)

    if not datasets:
        datasets = get_datasets(seqr_datasets_only=rd)

    excluded_sgs = set()
    if to_path(output_file).exists():
        logger.info(f'Output file {output_file} already exists. Reusing it.')
        excluded_sgs = read_existing_outputs_sgs(to_path(output_file))

    logger.info(f'Searching for SGs to add from datasets: {", ".join(datasets)}')

    # Get the sequencing groups for each dataset
    dataset_sgs = get_sequencing_groups_by_dataset(datasets, excluded_sgs)
    if dataset_sgs:
        logger.info(f'Found {len(dataset_sgs)} datasets with sequencing groups:')
        logger.info('Additional sequencing groups by dataset:')
        for dataset in dataset_sgs:
            logger.info(f'    {dataset} :: {len(dataset_sgs[dataset])} sequencing groups')
    else:
        logger.info('No additional sequencing groups found in the specified datasets.')
        sys.exit(0)

    # Get the samtools stats filepaths for each sequencing group
    samtools_stats_filepaths_by_sg = get_samtools_stats_filepaths(dataset_sgs)

    # Read the samtools stats files and get the statistics
    datasets_sg_stats = get_samtools_stats_for_sgs(dataset_sgs, samtools_stats_filepaths_by_sg)

    bucket_name = None
    if output_file.startswith('gs://'):
        bucket_name = to_path(output_file).parts[1]

    # Write the statistics to a file
    write_all_samtools_stats_to_file(
        datasets_sg_stats,
        bucket_name=bucket_name,
        output_file=output_file,
    )


if __name__ == '__main__':
    main()
