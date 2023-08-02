import re
import pytest
from pathlib import Path
from functools import cached_property

from cpg_workflows.batch import Batch
from cpg_workflows.jobs.picard import (
    get_intervals,
    markdup,
    vcf_qc,
    picard_collect_metrics,
    picard_hs_metrics,
    picard_wgs_metrics,
)

from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


class TestPicard:
    @cached_property
    def default_config(self) -> PipelineConfig:
        return PipelineConfig(
            workflow=WorkflowConfig(
                dataset='picard-test',
                access_level='test',
                sequencing_type='genome',
                check_inputs=False,
            ),
            images={
                'picard': 'picard_image:1.3.0',
            },
            other={
                'references': {
                    'broad': {
                        'ref_fasta': 'hg38_reference.fa',
                    }
                },
            },
        )

    def _setup(self, config: PipelineConfig, tmp_path: Path) -> Batch:
        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)

        return batch

    @pytest.mark.parametrize('scatter_count', [-1, 1, 10, 15, 20])
    @pytest.mark.parametrize(
        'source_intervals_path', ['source_intervals.txt', 'src_intrvls.txt']
    )
    def test_get_intervals(
        self, tmp_path: Path, scatter_count: int, source_intervals_path: str
    ):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test

        # Assert error for invalid scatter count
        if scatter_count <= 0:
            with pytest.raises(AssertionError, match=str(scatter_count)):
                job, _ = get_intervals(
                    b=batch,
                    scatter_count=scatter_count,
                    source_intervals_path=tmp_path / source_intervals_path,
                    job_attrs=None,
                    output_prefix=tmp_path / 'intervals',
                )
            return

        job, _ = get_intervals(
            b=batch,
            scatter_count=scatter_count,
            source_intervals_path=tmp_path / source_intervals_path,
            job_attrs=None,
            output_prefix=tmp_path / 'intervals',
        )
        cmd = get_command_str(job)

        # ---- Assertions
        if scatter_count == 1:
            assert job is None
            return

        assert re.search(f'INPUT=\S*{source_intervals_path}', cmd)

        for i in range(1, scatter_count + 1):
            assert re.search(f'temp_{i:04d}_of_{scatter_count}', cmd)
            assert re.search(f'/{i}.interval_list', cmd)

    def test_markdup(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        pass

        # ---- Assertions
        assert False

    def test_vcf_qc(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        pass

        # ---- Assertions
        assert False

    def test_picard_collect_metrics(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        pass

        # ---- Assertions
        assert False

    def test_picard_hs_metrics(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        pass

        # ---- Assertions
        assert False

    def test_picard_wgs_metrics(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        pass

        # ---- Assertions
        assert False

    # def test_creates_one_align_job(self, tmp_path: Path):
    #     # ---- Test setup
    #     config = default_config()
    #     set_config(config, tmp_path / 'config.toml')

    #     dataset_id = config.workflow.dataset
    #     batch = create_local_batch(tmp_path)
    #     sg = create_sequencing_group(
    #         dataset=dataset_id,
    #         sequencing_type=config.workflow.sequencing_type,
    #         alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
    #     )

    #     # ---- The job that we want to test
    #     _ = align(
    #         b=batch,
    #         sequencing_group=sg,
    #         extra_label=dataset_id,
    #         aligner=Aligner.DRAGMAP,
    #         markdup_tool=MarkDupTool.NO_MARKDUP,
    #     )

    #     # ---- Assertions
    #     align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
    #     assert len(align_jobs) == 1

    # def test_sorts_output_with_bamtools(self, tmp_path: Path):
    #     # ---- Test setup
    #     config = default_config()
    #     set_config(config.as_dict(), tmp_path / 'config.toml')

    #     dataset_id = config.workflow.dataset
    #     batch = create_local_batch(tmp_path)
    #     sg = create_sequencing_group(
    #         dataset=dataset_id,
    #         sequencing_type=config.workflow.sequencing_type,
    #         alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
    #     )

    #     # ---- The job that we want to test
    #     _ = align(
    #         b=batch,
    #         sequencing_group=sg,
    #         extra_label=dataset_id,
    #         aligner=Aligner.DRAGMAP,
    #         markdup_tool=MarkDupTool.NO_MARKDUP,
    #     )

    #     # ---- Assertions
    #     align_jobs = batch.select_jobs(rf'(.*){dataset_id}(.*)')
    #     cmd = get_command_str(align_jobs[0])
    #     assert re.search(r'\| samtools sort .* -Obam', cmd)
