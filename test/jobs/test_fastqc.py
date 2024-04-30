import re
from functools import cached_property
from pathlib import Path

import pytest

from cpg_utils.hail_batch import Batch
from cpg_workflows.filetypes import BamPath, FastqPath
from cpg_workflows.jobs.fastqc import fastqc

from .. import set_config
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from .helpers import get_command_str


class TestFastqc:
    @cached_property
    def default_config(self) -> PipelineConfig:
        return PipelineConfig(
            workflow=WorkflowConfig(
                dataset='fastq-test',
                access_level='test',
                sequencing_type='genome',
                check_inputs=False,
            ),
            images={'fastqc': 'fastqc_image:1.2.3'},
            references={
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                },
            },
            other={
                'resource_overrides': {},
            },
        )

    def _setup(self, config: PipelineConfig, tmp_path: Path) -> Batch:
        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)

        return batch

    @pytest.mark.parametrize('bam', ['input.bam', 'file.bam'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_fastqc_bam(self, tmp_path: Path, bam: str, job_attrs: dict):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        bam_path = BamPath(
            path=tmp_path / bam,
            index_path=tmp_path / f'{bam}.bai',
        )

        job = fastqc(
            b=batch,
            output_html_path=tmp_path / 'output.html',
            output_zip_path=tmp_path / 'output.zip',
            input_path=bam_path,
            subsample=True,
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert re.search(rf'{bam}', cmd)
        assert re.search(r'ln .*_fastqc.html \S+/\S+/out_html', cmd)
        assert re.search(r'ln .*_fastqc.zip \S+/\S+/out_zip', cmd)
        assert job_attrs.items() <= job.attributes.items()

    @pytest.mark.parametrize('fastq', ['input.fastq', 'file.fastq', 'in.fq.gz'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_fastqc_fastq(self, tmp_path: Path, fastq: str, job_attrs: dict):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        fastq_path: FastqPath = tmp_path / fastq

        job = fastqc(
            b=batch,
            output_html_path=tmp_path / 'output.html',
            output_zip_path=tmp_path / 'output.zip',
            input_path=fastq_path,
            subsample=True,
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert re.search(rf'{fastq}', cmd)
        assert re.search(r'ln .*_fastqc.html \S+/\S+/out_html', cmd)
        assert re.search(r'ln .*_fastqc.zip \S+/\S+/out_zip', cmd)
        assert job_attrs.items() <= job.attributes.items()
