"""
Test GATK-SV accessory functions.
"""

import pytest

from cpg_utils import to_path
from cpg_utils.config import set_config_paths
from cpg_workflows.stages.gatk_sv import (
    image_path,
    get_images,
    get_fasta,
    get_references,
)


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

    with open(tmp_path / 'config.toml', 'w') as fh:
        fh.write(TOML)
    set_config_paths([str(tmp_path / 'config.toml')])

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    image_dict = get_images(['first', 'second'])
    assert image_dict == {'first': 'first/image', 'second': 'second/image'}


def test_get_images_missing_allowed(tmp_path):
    """
    check that the image return works correctly
    """

    with open(tmp_path / 'config.toml', 'w') as fh:
        fh.write(TOML)
    set_config_paths([str(tmp_path / 'config.toml')])

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    image_dict = get_images(['first', 'second', 'third'], allow_missing=True)
    assert image_dict == {'first': 'first/image', 'second': 'second/image'}


def test_get_images_invalid(tmp_path):
    """
    check for failure when a requested key is missing
    """

    with open(tmp_path / 'config.toml', 'w') as fh:
        fh.write(TOML)
    set_config_paths([str(tmp_path / 'config.toml')])

    assert image_path('first') == 'first/image'
    assert image_path('second') == 'second/image'
    with pytest.raises(KeyError):
        get_images(['first', 'second', 'third'])


def test_get_fasta(tmp_path):
    """
    check that the fasta return works correctly
    """
    assert get_fasta() == to_path('this/is/the/ref.fasta')


def test_get_references(tmp_path):
    """
    check that the reference return works correctly
    """

    with open(tmp_path / 'config.toml', 'w') as fh:
        fh.write(TOML)
    set_config_paths([str(tmp_path / 'config.toml')])
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
