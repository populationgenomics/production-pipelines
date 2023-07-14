"""
concept test of a stage using a test_mode
"""
import inspect
import re
from unittest import mock

import pytest
from cpg_utils import to_path

from cpg_workflows.batch import get_batch
from cpg_workflows.stages.gatk_sv.gatk_sv_common import (
    _sv_batch_meta,
    _sv_individual_meta,
    add_gatk_sv_jobs,
    get_fasta,
    get_images,
    get_ref_panel,
    get_references,
    make_combined_ped
)

from test import set_config

from test.test_workflow import mock_create_cohort


TOML = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
driver_image = 'test'
sequencing_type = 'genome'
only_stages = ["{only}"]

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'
ref_fasta = 'strawberry'

[storage.default]
default = "{directory}"

[storage.fewgenomes]
default = "{directory}"

[storage.my_dataset]
default = "{directory}"

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'
dry_run = true

[references]
refgenome_sdf = "refgenome_sdf"
stratification = "{directory}"

[references.broad]
ref_fasta = "ref_fasta"

[images]
cpg_workflows = "cpgw_test"
"""


def test_sv_batch_meta():
    assert _sv_batch_meta('') == {'type': 'gatk-sv-batch-calls'}


def test_sv_individual_meta():
    assert _sv_individual_meta('') == {'type': 'gatk-sv-sequence-group-calls'}


def test_get_fasta(tmp_path):
    set_config(TOML, tmp_path / 'config.toml')
    assert get_fasta() == to_path('strawberry')


def test_get_images(tmp_path):

    set_config(TOML, tmp_path / 'config.toml')
    assert get_images(['cpg_workflows']) == {'cpg_workflows': 'cpgw_test'}


def test_get_images_allow_subset(tmp_path):

    set_config(TOML, tmp_path / 'config.toml')
    assert get_images(['cpg_workflows', 'piggly_wiggly'], True) == {'cpg_workflows': 'cpgw_test'}


def test_get_images_dont_allow_subset(tmp_path):

    set_config(TOML, tmp_path / 'config.toml')
    with pytest.raises(KeyError):
        get_images(['cpg_workflows', 'piggly_wiggly'], False)


@mock.patch('cpg_workflows.inputs.create_cohort', mock_create_cohort)
@mock.patch('cpg_workflows.workflow.list_of_all_dir_contents', set)
def test_sv_batch_meta(tmp_path):

    from cpg_workflows.workflow import run_workflow