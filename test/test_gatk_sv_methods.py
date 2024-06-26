"""
Test GATK-SV accessory functions.
"""

import pandas as pd
import pytest

from cpg_utils import to_path
from cpg_workflows.jobs.sample_batching import batch_sgs
from cpg_workflows.stages.gatk_sv.gatk_sv_common import get_fasta, get_images, get_references, image_path

from . import set_config

TOML = """
[workflow]
sequencing_type = 'genome'
ref_fasta = 'this/is/the/ref.fasta'

[images]
first = 'first/image'
second = 'second/image'

[references.gatk_sv]
sv_reference = 'gatk_sv_content'

[references.broad]
broad_reference = 'broad_content'
"""


def test_get_images(tmp_path):
    """
    check that the image return works correctly
    """
    set_config(TOML, tmp_path / 'config.toml')

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    image_dict = get_images(['first', 'second'])
    assert image_dict == {'first': 'first/image', 'second': 'second/image'}


def test_get_images_missing_allowed(tmp_path):
    """
    check that the image return works correctly
    """
    set_config(TOML, tmp_path / 'config.toml')

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    image_dict = get_images(['first', 'second', 'third'], allow_missing=True)
    assert image_dict == {'first': 'first/image', 'second': 'second/image'}


def test_get_images_invalid(tmp_path):
    """
    check for failure when a requested key is missing
    """
    set_config(TOML, tmp_path / 'config.toml')

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    with pytest.raises(KeyError):
        get_images(['first', 'second', 'third'])


def test_get_fasta(tmp_path):
    """
    check that the fasta return works correctly
    """
    set_config(TOML, tmp_path / 'config.toml')

    assert get_fasta() == to_path('this/is/the/ref.fasta')


def test_get_references(tmp_path):
    """
    check that the reference return works correctly
    """
    set_config(TOML, tmp_path / 'config.toml')

    assert get_references(['sv_reference', 'broad_reference']) == {
        'sv_reference': 'gatk_sv_content',
        'broad_reference': 'broad_content',
    }
    assert get_references([{'sneaky': 'sv_reference'}, 'broad_reference']) == {
        'sneaky': 'gatk_sv_content',
        'broad_reference': 'broad_content',
    }
    assert get_references([{'sneaky': 'sv_reference'}, 'weasel.broad_reference']) == {
        'sneaky': 'gatk_sv_content',
        'weasel.broad_reference': 'broad_content',
    }


def test_batch_samples():
    """
    use some dummy data to check that the batch_samples function works
    check for an appropriate size for each batch
    assert that each batch in turn has a higher mean coverage
    """

    # loading and formatting done in partition_batches function
    qc_table = to_path(__file__).parent / 'data' / 'gatk_sv' / 'evidence_qc.tsv'
    qc_df = pd.read_csv(qc_table, sep='\t')
    qc_df.columns = [x.replace('#', '') for x in qc_df.columns]

    # parameters to use
    min_size = 2
    max_size = 5

    # generate batches
    batches = batch_sgs(qc_df, min_size, max_size)

    # check that each incremental batch has higher mean coverage
    # and a size within range
    cov = 0
    for batch in batches:
        # test is max+1 as we're joining from male and female lists
        # separately, so sizes can be odd.
        # at ~hundreds of samples this is a fine margin of error
        assert min_size <= batch['size'] <= max_size + 1
        mean_of_medians = sum(batch['coverage_medians']) / len(batch['coverage_medians'])
        assert mean_of_medians > cov
        cov = mean_of_medians


def test_batch_samples_single_chunk():
    """
    use some dummy data to check that the batch_samples function works
    check for an appropriate size for each batch
    assert that each batch in turn has a higher mean coverage
    """

    # loading and formatting done in partition_batches function
    qc_table = to_path(__file__).parent / 'data' / 'gatk_sv' / 'evidence_qc.tsv'
    qc_df = pd.read_csv(qc_table, sep='\t')
    qc_df.columns = [x.replace('#', '') for x in qc_df.columns]

    # parameters to use
    min_size = 2
    max_size = 500

    # generate batches
    batches = batch_sgs(qc_df, min_size, max_size)
    assert len(batches) == 1
