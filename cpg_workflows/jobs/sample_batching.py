"""
extract from Broad sample batching script for GATK-SV
"""


import logging
import json

import numpy as np
import pandas as pd

from cpg_utils import to_path


SEX_VALS = {'male', 'female'}


def batch_samples(
    md: pd.DataFrame,
    batch_size: int = 200,
    min_batch_size: int = 100,
    max_batch_size: int = 300,
    cov_bins: int = 1,
):
    """
    Batch samples by coverage, dosage bias, and chrX ploidy

    Args:
        md (str): DataFrame of metadata
        batch_size (int): preferred batch size
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size
        cov_bins (int): number of coverage bins to use

    Returns:
        A dictionary
    """

    # Split samples based on >= 2 copies of chrX vs. < 2 copies
    isfemale = md.chrX_CopyNumber_rounded >= 2

    # Calculate approximate number of batches
    n_samples = len(md)

    # allow for shortcut?
    if min_batch_size <= n_samples <= max_batch_size:
        logging.info(f'Number of samples ({n_samples}) is within range of batch sizes')
        return {'1': md.ID.tolist()}

    elif n_samples < min_batch_size:
        raise ValueError(
            f'Number of samples ({n_samples}) is less than minimum batch size ({min_batch_size})'
        )

    n_per_cov = (n_samples // cov_bins) or 1
    batches_per_cov = (n_per_cov // batch_size) or 1
    n_per_batch = (n_per_cov // batches_per_cov) or 1

    logging.info(
        f"""
    n_samples: {n_samples}
    cov_bins: {cov_bins}
    n_per_cov: {n_per_cov}
    batches_per_cov: {batches_per_cov}
    n_per_batch: {n_per_batch}
"""
    )

    # Adjust batches to fit within specified range
    while n_per_batch < min_batch_size:
        batches_per_cov -= 1
        # stop trying to divide by 0
        if batches_per_cov == 0:
            raise ValueError('Batch size is too small for number of samples')
        n_per_batch = n_per_cov // batches_per_cov

    while n_per_batch > max_batch_size:
        batches_per_cov += 1
        # stop trying to divide by 0
        if batches_per_cov == 0:
            raise ValueError('Batch size is too small for number of samples')
        n_per_batch = n_per_cov // batches_per_cov

    # Split samples by sex, then roughly by coverage
    md_sex = {
        'male': md[~isfemale].sort_values(by='median_coverage'),
        'female': md[isfemale].sort_values(by='median_coverage'),
    }
    md_sex_cov = {sex: np.array_split(md_sex[sex], cov_bins) for sex in SEX_VALS}

    # Create batches & assign sample IDs
    i_vals = range(batches_per_cov)
    batches = {}
    for cov in range(cov_bins):
        # Order by WGD within each subset
        submds = {
            sex: np.array_split(
                md_sex_cov[sex][cov].sort_values(by='wgd_score'), batches_per_cov
            )
            for sex in SEX_VALS
        }
        submds = [pd.concat([submds['male'][i], submds['female'][i]]) for i in i_vals]
        for i in i_vals:
            batches[f'c{cov}b{i}'] = submds[i].ID.tolist()

    return batches


def partition_batches(metadata_file: str, output_json: str):
    """
    Runs this process
    - loads the input data
    - obtains batches
    - writes out the data into one single cross-batch file
    - also writes out per-batch sample lists

    Args:
        metadata_file (str): path to the metadata file
        output_json (str): location to write the batched samples outå
    """

    # read in the metadata contents
    logging.basicConfig(level=logging.INFO)

    # load in the metadata file
    md = pd.read_csv(metadata_file, sep='\t', low_memory=False)
    md.columns = [x.replace('#', '') for x in md.columns]
    batches = batch_samples(md)

    # write out the batches to GCP
    logging.info(batches)

    with to_path(output_json).open('w') as handle:
        json.dump(batches, handle, indent=4)
