"""
extract from Broad sample batching script for GATK-SV
"""


import logging
import json

import numpy as np
import pandas as pd

from cpg_utils import to_path


SEX_VALS = {'male', 'female'}
PCR_VALS = {'PCRPLUS', 'PCRMINUS'}


def batch_samples(md: pd.DataFrame, min_batch_size, max_batch_size) -> list[dict]:
    """
    Batch samples by coverage, dosage bias, and chrX ploidy

    Args:
        md (str): DataFrame of metadata
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size

    Returns:
        A list of batches, each with samples to fit into each batch
        {
            ID: {
                'samples': [sample1, sample2, ...],
                'male/female': float,
                'size': int,
            }
        }
    """

    # Split samples based on >= 2 copies of chrX vs. < 2 copies
    isfemale = md.chrX_CopyNumber_rounded >= 2
    ismale = md.chrX_CopyNumber_rounded == 1

    # Calculate approximate number of batches
    n_samples = len(md)

    # allow for shortcut?
    if min_batch_size <= n_samples <= max_batch_size:
        logging.info(f'Number of samples ({n_samples}) is within range of batch sizes')
        return [
            {
                'samples': md.ID.tolist(),
                'mf_ratio': ismale / isfemale,
                'size': n_samples,
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
        'male': md[~isfemale].sort_values(by='median_coverage'),
        'female': md[isfemale].sort_values(by='median_coverage'),
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
            }
        )

    return batches


def partition_batches(
    metadata_file: str,
    output_json: str,
    batch_size: int,
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
        output_json (str): location to write the batched samples out
        batch_size (int): preferred batch size
        min_batch_size (int): minimum batch size
        max_batch_size (int): maximum batch size
    """

    # read in the metadata contents
    logging.basicConfig(level=logging.INFO)

    # load in the metadata file
    md = pd.read_csv(metadata_file, sep='\t', low_memory=False)
    md.columns = [x.replace('#', '') for x in md.columns]

    # split 'em up into inferred PCR +/-
    for pcr in PCR_VALS:
        pcr_state_samples = md[(md.inferred_pcr_status == pcr)]
        if len(pcr_state_samples) < min_batch_size:
            logging.info(f'Insufficient samples found for PCR status {pcr}')
            continue

        batches = batch_samples(pcr_state_samples, min_batch_size, max_batch_size)

        # write out the batches to GCP
        logging.info(batches)

        # write the batch JSON out to a file, inc. PCR state
        with to_path(output_json.format(pcr)).open('w') as handle:
            json.dump(batches, handle, indent=4)
