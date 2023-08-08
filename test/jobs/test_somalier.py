import pandas as pd
from cpg_utils import Path
from hailtop.batch import Batch, Resource
from hailtop.batch.job import Job
from pytest_mock import MockFixture

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='som-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'somalier': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                    'somalier_sites': '',
                }
            },
            'cramqc': {'num_pcs': '4'},
        },
    )
