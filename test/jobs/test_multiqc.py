import os
import re
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.jobs.samtools import samtools_stats

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='samtools-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'samtools': 'test_image',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
            }
        },
    )


def setup_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    batch = create_local_batch(tmp_path)
    cram_pth = create_cram_input(location=tmp_path, index=True)

    return config, cram_pth, batch
