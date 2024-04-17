import random
import re
from functools import cached_property
from pathlib import Path

import pytest

from cpg_utils.hail_batch import Batch, fasta_res_group
from cpg_workflows.filetypes import BamPath, CramPath, GvcfPath
from cpg_workflows.jobs.picard import (
    get_intervals,
    markdup,
    picard_collect_metrics,
    picard_hs_metrics,
    picard_wgs_metrics,
    vcf_qc,
)

from .. import set_config
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from .helpers import get_command_str, get_path_from_resource_file


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
            references={
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dbsnp_vcf': 'dbsnp.vcf.gz',
                    'dbsnp_vcf_index': 'dbsnp.vcf.gz.tbi',
                    'genome_evaluation_interval_lists': 'intervals.txt',
                    'exome_evaluation_interval_lists': 'exome_intervals.txt',
                },
                'hg38_telomeres_and_centromeres_intervals': {
                    'interval_list': 'hg38_telomeres_and_centromeres.interval_list',
                },
            },
            other={
                'resource_overrides': {
                    'picard_storage_gb': 1,
                },
                'cramqc': {'assume_sorted': False},
            },
        )

    def _setup(self, config: PipelineConfig, tmp_path: Path, ref: str | None = None) -> Batch:
        if ref is not None:
            config.workflow.ref_fasta = ref

        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)

        return batch

    @pytest.mark.parametrize('scatter_count', [-645, -1, 0])
    def test_get_intervals_scatter_less_than_zero(self, tmp_path: Path, scatter_count: int):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test

        # Assert error for invalid scatter count
        with pytest.raises(AssertionError, match=str(scatter_count)):
            _, _ = get_intervals(
                b=batch,
                scatter_count=scatter_count,
                source_intervals_path=tmp_path / 'source_intervals.txt',
                output_prefix=tmp_path / 'intervals',
            )

    def test_get_intervals_scatter_is_one(self, tmp_path: Path):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        job, _ = get_intervals(
            b=batch,
            scatter_count=1,
            source_intervals_path=tmp_path / 'source_intervals.txt',
            output_prefix=tmp_path,
        )

        # ---- Assertions
        assert job is None

    @pytest.mark.parametrize('scatter_count', [10, 15, 20])
    def test_get_intervals_all_outputs_exist(
        self,
        tmp_path: Path,
        scatter_count: int,
    ):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # Create all of the files
        output_prefix = tmp_path
        for i in range(1, scatter_count + 1):
            file = output_prefix / f'{i}.interval_list'
            file.touch()

        # ---- The job we want to test
        job, intervals = get_intervals(
            b=batch,
            scatter_count=scatter_count,
            source_intervals_path=tmp_path / 'source_intervals.txt',
            output_prefix=output_prefix,
        )

        # ---- Assertions
        assert job is None
        assert len(intervals) == scatter_count

        for i, file in enumerate(intervals):
            file_path = get_path_from_resource_file(file)
            assert file_path == str(tmp_path / f'{i+1}.interval_list')

    @pytest.mark.parametrize('scatter_count', [10, 15, 20])
    def test_get_intervals_some_outputs_exist(
        self,
        tmp_path: Path,
        scatter_count: int,
    ):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # Create only some of the files
        output_prefix = tmp_path
        for i in range(1, random.randint(2, scatter_count)):
            file = output_prefix / f'{i}.interval_list'
            file.touch()

        # ---- The job we want to test
        job, intervals = get_intervals(
            b=batch,
            scatter_count=scatter_count,
            source_intervals_path=tmp_path / 'source_intervals.txt',
            output_prefix=output_prefix,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job is not None
        assert len(intervals) == scatter_count

        for i in range(1, scatter_count + 1):
            assert re.search(rf'temp_{i:04d}_of_{scatter_count}', cmd)
            assert re.search(rf'/{i}.interval_list', cmd)

    @pytest.mark.parametrize('scatter_count', [10, 15, 20])
    @pytest.mark.parametrize('source_intervals_path', ['source_intervals.txt', 'src_intrvls.txt'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_get_intervals_valid_inputs(
        self,
        tmp_path: Path,
        scatter_count: int,
        source_intervals_path: str,
        job_attrs: dict,
    ):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        job, _ = get_intervals(
            b=batch,
            scatter_count=scatter_count,
            source_intervals_path=tmp_path / source_intervals_path,
            output_prefix=tmp_path / 'intervals',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'-I \S*{source_intervals_path}', cmd)

        for i in range(1, scatter_count + 1):
            assert re.search(rf'temp_{i:04d}_of_{scatter_count}', cmd)
            assert re.search(rf'/{i}.interval_list', cmd)

    # test_get_intervals_already_exists
    # put files into data folder

    @pytest.mark.parametrize('ref', [None, 'ref.fa'])
    def test_markdup_ref_fasta(self, tmp_path: Path, ref: str):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path, ref=ref)

        # ---- The job we want to test
        sorted_bam = BamPath(path=tmp_path / 'in.bam', index_path=tmp_path / 'in.bam.bai')
        job = markdup(
            b=batch,
            sorted_bam=sorted_bam,
            fasta_reference=fasta_res_group(batch),
        )
        cmd = get_command_str(job)

        # ---- Assertions
        if ref is None:
            ref = self.default_config.references['broad']['ref_fasta']

        assert re.search(rf'-T\s*\S*{ref}', cmd)

    @pytest.mark.parametrize('bam', ['input.bam', 'sorted.bam'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_markdup_valid_inputs(self, tmp_path: Path, bam: str, job_attrs: dict):
        # ---- Test setup

        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        sorted_bam = BamPath(path=tmp_path / bam, index_path=tmp_path / f'{bam}.bai')
        job = markdup(
            b=batch,
            sorted_bam=sorted_bam,
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'I=\S*{bam}', cmd)

    @pytest.mark.parametrize('gvcf', ['file.gvcf.gz', 'file.gvcf'])
    @pytest.mark.parametrize('dbsnp', ['dbsnp.vcf.gz', 'DBSNP.vcf.gz'])
    @pytest.mark.parametrize('intervals', ['intervals.txt', 'intrvls.txt'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_vcf_qc(self, tmp_path: Path, gvcf: str, dbsnp: str, intervals: str, job_attrs: dict):
        # ---- Test setup
        config = self.default_config
        config.references['broad'] = {
            'ref_fasta': 'hg38_reference.fa',
            'dbsnp_vcf': dbsnp,
            'dbsnp_vcf_index': dbsnp + '.tbi',
            'genome_evaluation_interval_lists': intervals,
        }
        batch = self._setup(config, tmp_path)

        # ---- The job we want to test
        gvcf_path = GvcfPath(tmp_path / gvcf)
        job = vcf_qc(
            b=batch,
            vcf_or_gvcf=gvcf_path.resource_group(batch),
            is_gvcf=True,
            sequencing_group_count=100,
            output_summary_path=tmp_path / 'summary.txt',
            output_detail_path=tmp_path / 'output.txt',
            overwrite=False,
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'INPUT=\S*{gvcf}', cmd)
        assert re.search(rf'DBSNP=\S*{dbsnp}', cmd)
        assert re.search(rf'TARGET_INTERVALS=\S*{intervals}', cmd)
        assert re.search(r'GVCF_INPUT=true', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize('assume_sorted', [True, False])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_collect_metrics(self, tmp_path: Path, cram: str, assume_sorted: bool, job_attrs: dict):
        # ---- Test setup
        config = self.default_config
        config.other['cramqc']['assume_sorted'] = assume_sorted
        batch = self._setup(config, tmp_path)

        # ---- The job we want to test
        cram_path = CramPath(path=tmp_path / cram, index_path=tmp_path / (cram + '.crai'))
        job = picard_collect_metrics(
            b=batch,
            cram_path=cram_path,
            out_alignment_summary_metrics_path=tmp_path / 'algn_summary_metrics.txt',
            out_base_distribution_by_cycle_metrics_path=tmp_path / 'base_dist.txt',
            out_insert_size_metrics_path=tmp_path / 'insrt_size_metrics.txt',
            out_quality_by_cycle_metrics_path=tmp_path / 'qual_by_cycle.txt',
            out_quality_yield_metrics_path=tmp_path / 'qual_yield_by_cycle.txt',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'ASSUME_SORTED={assume_sorted}', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize('exome_intervals', ['exome_intervals.txt', 'ex_intrvls.txt'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_hs_metrics(self, tmp_path: Path, cram: str, exome_intervals: str, job_attrs: dict):
        # ---- Test setup
        config = self.default_config
        config.workflow = WorkflowConfig(
            dataset='picard-test',
            access_level='test',
            sequencing_type='exome',
            check_inputs=False,
        )
        config.references['broad'] = {
            'ref_fasta': 'hg38_reference.fa',
            'exome_evaluation_interval_lists': exome_intervals,
        }
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        cram_path = CramPath(path=tmp_path / cram, index_path=tmp_path / (cram + '.crai'))
        job = picard_hs_metrics(
            b=batch,
            cram_path=cram_path,
            out_picard_hs_metrics_path=tmp_path / 'picard_hs_metrics.txt',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'I=\S+{exome_intervals}', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize('genome_intervals', ['genome_intervals.txt', 'gnm_intrvls.txt'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_wgs_metrics(self, tmp_path: Path, cram: str, genome_intervals: str, job_attrs: dict):
        # ---- Test setup
        config = self.default_config
        config.references['broad'] = {
            'ref_fasta': 'hg38_reference.fa',
            'genome_coverage_interval_list': genome_intervals,
        }
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        cram_path = CramPath(path=tmp_path / cram, index_path=tmp_path / (cram + '.crai'))
        job = picard_wgs_metrics(
            b=batch,
            cram_path=cram_path,
            out_picard_wgs_metrics_path=tmp_path / 'picard_hs_metrics.txt',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'INTERVALS=\S+{genome_intervals}', cmd)
