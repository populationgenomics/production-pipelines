import re

import pytest
from cpg_utils import Path
from pytest_mock import MockFixture

from cpg_workflows.jobs.multiqc import (
    _write_sg_id_map,
    check_multiqc,
    check_report_job,
    multiqc,
)

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.dataset import create_dataset
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='multiqc-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
            dry_run=True,  # may need to change this
        ),
        images={
            'multiqc': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                },
            },
        },
    )


def setup_multiqc_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    dataset_id = config.workflow.dataset

    dataset = create_dataset(name=dataset_id)
    dataset.add_sequencing_group(id='CPG000001', external_id='SAMPLE1')
    batch = create_local_batch(tmp_path)

    paths = [(tmp_path / f'{i}') for i in range(4)]
    return config, batch, dataset, paths


class TestMultiQC:
    def test_creates_one_multiqc_job(self, tmp_path: Path):
        # ---- Test setup
        config, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )

        # ---- Assertions
        multiqc_job = batch.select_jobs('MultiQC')
        assert len(multiqc_job) == 1
