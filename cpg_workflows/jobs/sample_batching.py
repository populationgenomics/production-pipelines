"""
extract from Broad sample batching script for GATK-SV
"""


import logging
import json

import numpy as np
import pandas as pd

from cpg_utils import to_path


SEX_VALS = {'male', 'female'}


def batch_samples(md: pd.DataFrame, min_batch_size, max_batch_size) -> list[dict]:
    """
    Batch samples by coverage, and chrX ploidy

    Starts with a base assumption of one batch, increasing the batch count
    until min < batch size < max

    When an appropriate number of batches is found, the samples are split
    by sex (assigned by inferred X ploidy), then all samples are ordered
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
        A list of batches, each with samples
        Each batch self-documents its size and male/female ratio
        {
            batch_ID: {
                'samples': [sample1, sample2, ...],
                'mf_ratio': float,
                'size': int,
                'coverage_medians': [float, float, ...],
            }
        }
    """

    # Split samples based on >= 2 copies of chrX vs. < 2 copies
    is_female = md.chrX_CopyNumber_rounded >= 2
    is_male = md.chrX_CopyNumber_rounded == 1

    # region: Calculate approximate number of batches
    n_samples = len(md)

    # shortcut return if the total samples is a valid batch
    if min_batch_size <= n_samples <= max_batch_size:
        logging.info(f'Number of samples ({n_samples}) is within range of batch sizes')
        return [
            {
                'samples': md.ID.tolist(),
                'mf_ratio': is_male / is_female,
                'size': n_samples,
                'coverage_medians': md.median_coverage.tolist(),
            }
        ]

    # start with a single bin and work upwards
    cov_bins = 1
    n_per_batch = n_samples // cov_bins

    # add more bins until we get to a reasonable number of samples per batch
    # shrink in one direction by increasing num. batches
    while n_per_batch > max_batch_size:
        cov_bins += 1
        n_per_batch = n_samples // cov_bins
    # endregion

    logging.info(
        f"""
    Batching Specs:
    Total samples: {n_samples}
    Num. bins: {cov_bins}
    Approx samples per batch: {n_per_batch}
"""
    )

    # Split samples by sex, then roughly by coverage
    md_sex = {
        'male': md[~is_female].sort_values(by='median_coverage'),
        'female': md[is_female].sort_values(by='median_coverage'),
    }

    # array split each sex proportionally
    md_sex_cov = {sex: np.array_split(md_sex[sex], cov_bins) for sex in SEX_VALS}

    # create batches
    batches = []
    for cov in range(cov_bins):
        sample_ids = pd.concat([md_sex_cov['male'][cov], md_sex_cov['female'][cov]])
        batches.append(
            {
                'samples': sample_ids.ID.tolist(),
                'size': len(sample_ids),
                'mf_ratio': len(md_sex_cov['male'][cov])
                / len(md_sex_cov['female'][cov]),
                'coverage_medians': sample_ids.median_coverage.tolist(),
            }
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
        output_json (str): location to write the batched samples out
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size
    """

    # read in the metadata contents
    logging.basicConfig(level=logging.INFO)

    # load in the metadata file
    md = pd.read_csv(metadata_file, sep='\t', low_memory=False)
    md.columns = [x.replace('#', '') for x in md.columns]

    # filter to the PCR-state samples we're interested in
    # surely there's a neater way to do this...
    md = md[np.array([sam in sample_ids for sam in md.ID.tolist()])]

    # check that we have enough samples to batch
    # should have already been checked prior to Stage starting
    if len(md) < min_batch_size:
        raise ValueError('Insufficient samples found for batch generation')

    # generate the batches
    batches = batch_samples(
        md=md, min_batch_size=min_batch_size, max_batch_size=max_batch_size
    )

    # write out the batches to GCP
    logging.info(batches)

    # write the batch JSON out to a file, inc. PCR state
    with to_path(output_json).open('w') as handle:
        json.dump(batches, handle, indent=4)
