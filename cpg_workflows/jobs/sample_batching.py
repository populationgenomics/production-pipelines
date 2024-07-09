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
                'sequencing_groups': [sample1, sample2, ...],
                'male': [sample1, sample3, ...],
                'female': [sample1, sample3, ...],
                'mf_ratio': float,
                'size': int,
                'coverage_medians': [float, float, ...],
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
                'sequencing_groups': md.ID.tolist(),
                'male': md[~is_female].ID.to_list(),
                'male_count': len(md[~is_female]),
                'female': md[is_female].ID.to_list(),
                'female_count': len(md[is_female]),
                'mf_ratio': is_male.sum() / is_female.sum(),
                'size': n_sg,
                'coverage_medians': md.median_coverage.tolist(),
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
                'sequencing_groups': sample_ids.ID.tolist(),
                'size': len(sample_ids),
                'male': md_sex_cov['male'][cov].ID.tolist(),  # type: ignore
                'male_count': len(md_sex_cov['male'][cov]),
                'female': md_sex_cov['female'][cov].ID.tolist(),  # type: ignore
                'female_count': len(md_sex_cov['female'][cov]),
                'mf_ratio': (len(md_sex_cov['male'][cov]) or 1) / (len(md_sex_cov['female'][cov]) or 1),
                'coverage_medians': sample_ids.median_coverage.tolist(),
            },
        )

    return batches


def partition_batches(
    metadata_file: str,
    sample_ids: list[str],
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
        metadata_file (str): path to the metadata file
        sample_ids (list[str]): sample IDs to consider
        output_json (str): location to write the batch result
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size
    """

    # read in the metadata contents
    get_logger(__file__).info('Starting the batch creation process')

    # load in the metadata file
    md = pd.read_csv(metadata_file, sep='\t', low_memory=False)
    md.columns = [x.replace('#', '') for x in md.columns]  # type: ignore

    # filter to the PCR-state SGs we're interested in
    md = md.query('ID in @sample_ids')

    # check that we have enough samples to batch
    # should have already been checked prior to Stage starting
    if len(md) < min_batch_size:
        raise ValueError('Insufficient Seq Groups found for batch generation')

    # generate the batches
    batches = batch_sgs(md=md, min_batch_size=min_batch_size, max_batch_size=max_batch_size)

    # write out the batches to GCP
    get_logger().info(batches)

    # write the batch JSON out to a file, inc. PCR state
    with to_path(output_json).open('w') as handle:
        json.dump(batches, handle, indent=4)
