"""
extract from Broad sample batching script for GATK-SV
"""

import json

import numpy as np
import pandas as pd

from cpg_utils import to_path
from cpg_workflows.utils import get_logger

SEX_VALS = {'male', 'female'}


def batch_sgs(md: pd.DataFrame, min_batch_size: int, max_batch_size: int) -> list[dict]:
    """
    Batch sequencing groups by coverage, and chrX ploidy

    Starts with a base assumption of one batch, increasing the batch count
    until min < batch size < max

    When an appropriate number of batches is found, the SGs are split
    by sex (assigned by inferred X ploidy), then all SGs are ordered
    by median coverage.

    Using numpy.array_split, the ordered male & female sample lists are split
    into equal sized sub-lists, then the sub-lists are combined into batches.
    This groups lowest coverage male and female samples into a single batch,
    continuing up to the highest coverage batch.

    This method maintains a consistent sex ratio across all batches

    Args:
        md (str): DataFrame of metadata
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size

    Returns:
        A list of batches, each with sequencing_groups
        Each batch self-documents its size and male/female ratio
        {
            batch_ID: {
                'size': int,
                'male_count': int,
                'female_count': int,
                'mf_ratio': float,
                'batch_median_coverage': float,
                'coverage_range': (float, float),
                'sequencing_groups': [sample1, sample2, ...],
                'coverage_medians': [float, float, ...],
                'male': [sample1, sample3, ...],
                'female': [sample1, sample3, ...],
            }
        }
    """

    # Split SGs based on >= 2 copies of chrX vs. < 2 copies
    is_female = md.chrX_CopyNumber_rounded >= 2
    is_male = md.chrX_CopyNumber_rounded == 1

    # region: Calculate approximate number of batches
    n_sg = len(md)

    # shortcut return if the total sequencing groups is a valid batch
    if min_batch_size <= n_sg <= max_batch_size:
        get_logger().info(f'Number of sequencing_groups ({n_sg}) is within range of batch sizes')
        return [
            {
                'size': n_sg,
                'male_count': len(md[~is_female]),
                'female_count': len(md[is_female]),
                'mf_ratio': is_male.sum() / is_female.sum(),
                'batch_median_coverage': md.median_coverage.median(),
                'coverage_range': (md.median_coverage.min(), md.median_coverage.max()),
                'sequencing_groups': md.ID.tolist(),
                'coverage_medians': md.median_coverage.tolist(),
                'male': md[~is_female].ID.to_list(),
                'female': md[is_female].ID.to_list(),
            },
        ]

    # start with a single bin and work upwards
    cov_bins = 1
    n_per_batch = n_sg // cov_bins

    # add more bins until we get to a reasonable number of SGs per batch
    # shrink in one direction by increasing num. batches
    while n_per_batch > max_batch_size:
        cov_bins += 1
        n_per_batch = n_sg // cov_bins
    # endregion

    get_logger().info(
        f"""
    Batching Specs:
    Total sequencing_groups: {n_sg}
    Num. bins: {cov_bins}
    Approx sequencing groups per batch: {n_per_batch}
""",
    )

    # Split sequencing groups by sex, then roughly by coverage
    md_sex = {
        'male': md[~is_female].sort_values(by='median_coverage'),
        'female': md[is_female].sort_values(by='median_coverage'),
    }

    get_logger().info('Sex distribution across all samples:')
    for sex, entries in md_sex.items():
        get_logger().info(f'{sex}: {len(entries)}')

    # array split each sex proportionally
    md_sex_cov = {sex: np.array_split(md_sex[sex], cov_bins) for sex in SEX_VALS}

    # create batches
    batches = []
    for cov in range(cov_bins):
        sample_ids = pd.concat([md_sex_cov['male'][cov], md_sex_cov['female'][cov]])  # type: ignore
        get_logger().info(
            f"""
        ---
        Batch {cov}:
        Samples: {len(sample_ids)}
        Males: {len(md_sex_cov['male'][cov])}
        Females: {len(md_sex_cov['female'][cov])}
        ---
        """,
        )
        batches.append(
            {
                'size': len(sample_ids),
                'male_count': len(md_sex_cov['male'][cov]),
                'female_count': len(md_sex_cov['female'][cov]),
                'mf_ratio': (len(md_sex_cov['male'][cov]) or 1) / (len(md_sex_cov['female'][cov]) or 1),
                'batch_median_coverage': sample_ids.median_coverage.median(),
                'coverage_range': (sample_ids.median_coverage.min(), sample_ids.median_coverage.max()),
                'sequencing_groups': sample_ids.ID.tolist(),
                'coverage_medians': sample_ids.median_coverage.tolist(),
                'male': md_sex_cov['male'][cov].ID.tolist(),  # type: ignore
                'female': md_sex_cov['female'][cov].ID.tolist(),  # type: ignore
            },
        )

    return batches


def batch_sgs_by_library(md: pd.DataFrame, min_batch_size: int, max_batch_size: int) -> list[dict]:
    """
    Batch sequencing groups by coverage, chrX ploidy, and sequencing library

    Args:
        md (pd.DataFrame): DataFrame of metadata
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size

    Returns:
        A list of batches, each with sequencing_groups and additional information
    """
    # Check that the library column exists
    if 'library' not in md.columns:
        return batch_sgs(md, min_batch_size, max_batch_size)

    # Group SGs by library type
    library_groups = md.groupby('library')

    # Create batches from each library's DataFrame
    batches = []
    for _, library_df in library_groups:
        batches.extend(batch_sgs(library_df, min_batch_size, max_batch_size))

    return batches


def add_sg_meta_fields(sg_df: pd.DataFrame, sg_meta: dict[str, dict]) -> pd.DataFrame:
    """
    Adds the sg meta fields 'library' and 'facility' to the DataFrame

    Args:
        sg_df (pd.DataFrame): DataFrame of sequencing groups
        sg_meta (dict[str, dict]): meta dict keyed by SG ID

    Returns:
        pd.DataFrame: DataFrame with the sg_meta fields parsed and added
    """
    sg_df['library'] = sg_df['ID'].map(
        lambda x: sg_meta[x].get('library_type', sg_meta[x].get('sequencing_library', 'unknown')),
    )
    sg_df['facility'] = sg_df['ID'].map(
        lambda x: sg_meta[x].get('facility', sg_meta[x].get('sequencing_facility', 'unknown')),
    )
    return sg_df


def partition_batches(
    metadata_files: list[str],
    sequencing_groups_json: str,
    output_json: str,
    min_batch_size: int,
    max_batch_size: int,
):
    """
    Runs this process
    - loads the input data
    - obtains batches
    - writes out the data into one single cross-batch file
    - also writes out per-batch sample lists

    Args:
        metadata_files (list[str]): paths to the metadata files
        sequencing_groups_json (str): path to json of all SGs with their meta
        output_json (str): location to write the batch result
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size
    """

    # read in the metadata contents
    get_logger(__file__).info('Starting the batch creation process')

    # load in the metadata files
    md = pd.concat([pd.read_csv(md_file, sep='\t', low_memory=False) for md_file in metadata_files])
    md.columns = [x.replace('#', '') for x in md.columns]  # type: ignore

    # load the sequencing groups
    with to_path(sequencing_groups_json).open('r') as f:
        sequencing_groups = json.load(f)

    # filter to the PCR-state SGs we're interested in
    sample_ids = list(sequencing_groups.keys())
    md = md.query('ID in @sample_ids')
    md = add_sg_meta_fields(md, sequencing_groups)

    # check that we have enough samples to batch
    # should have already been checked prior to Stage starting
    if len(md) < min_batch_size:
        raise ValueError('Insufficient Seq Groups found for batch generation')

    # generate the batches
    batches = batch_sgs_by_library(md=md, min_batch_size=min_batch_size, max_batch_size=max_batch_size)

    # write out the batches to GCP
    get_logger().info(batches)

    # write the batch JSON out to a file, inc. PCR state
    with to_path(output_json).open('w') as handle:
        json.dump(batches, handle, indent=4)
