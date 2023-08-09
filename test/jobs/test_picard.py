import re
import pytest
from pathlib import Path
from functools import cached_property

from cpg_workflows.batch import Batch
from cpg_workflows.filetypes import BamPath, CramPath, GvcfPath
from cpg_workflows.jobs.picard import (
    get_intervals,
    markdup,
    vcf_qc,
    picard_collect_metrics,
    picard_hs_metrics,
    picard_wgs_metrics,
)

from .. import set_config
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
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
            references={
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dbsnp_vcf': 'dbsnp.vcf.gz',
                    'dbsnp_vcf_index': 'dbsnp.vcf.gz.tbi',
                    'genome_evaluation_interval_lists': 'intervals.txt',
                    'exome_evaluation_interval_lists': 'exome_intervals.txt',
                }
            },
            other={
                'resource_overrides': {
                    'picard_storage_gb': 1,
                },
                'cramqc': {'assume_sorted': False},
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
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_get_intervals(
        self,
        tmp_path: Path,
        scatter_count: int,
        source_intervals_path: str,
        job_attrs: dict,
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
                    output_prefix=tmp_path / 'intervals',
                    job_attrs=job_attrs,
                )
                assert job_attrs.items() <= job.attributes.items()
            return

        job, _ = get_intervals(
            b=batch,
            scatter_count=scatter_count,
            source_intervals_path=tmp_path / source_intervals_path,
            output_prefix=tmp_path / 'intervals',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        if scatter_count == 1:
            assert job is None
            return

        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'INPUT=\S*{source_intervals_path}', cmd)

        for i in range(1, scatter_count + 1):
            assert re.search(rf'temp_{i:04d}_of_{scatter_count}', cmd)
            assert re.search(rf'/{i}.interval_list', cmd)

    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_markdup(self, tmp_path: Path, job_attrs: dict):
        # ---- Test setup
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        sorted_bam = BamPath(
            path=tmp_path / 'sorted.bam', index_path=tmp_path / 'sorted.bam.bai'
        )
        job = markdup(
            b=batch,
            sorted_bam=sorted_bam,
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(r'I=\S*sorted.bam', cmd)
        assert re.search(r'-T\s*\S*hg38_reference.fa', cmd)

    @pytest.mark.parametrize('gvcf', ['file.gvcf.gz', 'file.gvcf'])
    @pytest.mark.parametrize('dbsnp', ['dbsnp.vcf.gz', 'DBSNP.vcf.gz'])
    @pytest.mark.parametrize('intervals', ['intervals.txt', 'intrvls.txt'])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_vcf_qc(
        self, tmp_path: Path, gvcf: str, dbsnp: str, intervals: str, job_attrs: dict
    ):
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
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'INPUT=\S*{gvcf}', cmd)
        assert re.search(rf'DBSNP=\S*{dbsnp}', cmd)
        assert re.search(rf'TARGET_INTERVALS=\S*{intervals}', cmd)
        assert re.search(rf'GVCF_INPUT=true', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize('assume_sorted', [True, False])
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_collect_metrics(
        self, tmp_path: Path, cram: str, assume_sorted: bool, job_attrs: dict
    ):
        # ---- Test setup
        config = self.default_config
        config.other['cramqc']['assume_sorted'] = assume_sorted
        batch = self._setup(config, tmp_path)

        # ---- The job we want to test
        cram_path = CramPath(
            path=tmp_path / cram, index_path=tmp_path / (cram + '.crai')
        )
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
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'ASSUME_SORTED={assume_sorted}', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize(
        'exome_intervals', ['exome_intervals.txt', 'ex_intrvls.txt']
    )
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_hs_metrics(
        self, tmp_path: Path, cram: str, exome_intervals: str, job_attrs: dict
    ):
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
        cram_path = CramPath(
            path=tmp_path / cram, index_path=tmp_path / (cram + '.crai')
        )
        job = picard_hs_metrics(
            b=batch,
            cram_path=cram_path,
            out_picard_hs_metrics_path=tmp_path / 'picard_hs_metrics.txt',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'I=\S+{exome_intervals}', cmd)

    @pytest.mark.parametrize('cram', ['file.cram', 'file2.cram'])
    @pytest.mark.parametrize(
        'genome_intervals', ['genome_intervals.txt', 'gnm_intrvls.txt']
    )
    @pytest.mark.parametrize('job_attrs', [{'blah': 'abc'}, {'test': '123'}])
    def test_picard_wgs_metrics(
        self, tmp_path: Path, cram: str, genome_intervals: str, job_attrs: dict
    ):
        # ---- Test setup
        config = self.default_config
        config.references['broad'] = {
            'ref_fasta': 'hg38_reference.fa',
            'genome_coverage_interval_list': genome_intervals,
        }
        batch = self._setup(self.default_config, tmp_path)

        # ---- The job we want to test
        cram_path = CramPath(
            path=tmp_path / cram, index_path=tmp_path / (cram + '.crai')
        )
        job = picard_wgs_metrics(
            b=batch,
            cram_path=cram_path,
            out_picard_wgs_metrics_path=tmp_path / 'picard_hs_metrics.txt',
            job_attrs=job_attrs,
        )
        cmd = get_command_str(job)

        # ---- Assertions
        assert job_attrs.items() <= job.attributes.items()
        assert re.search(rf'CRAM=\$BATCH_TMPDIR/{cram}', cmd)
        assert re.search(rf'CRAI=\$BATCH_TMPDIR/{cram}.crai', cmd)
        assert re.search(r'REFERENCE_SEQUENCE=\S+hg38_reference.fa', cmd)
        assert re.search(rf'INTERVALS=\S+{genome_intervals}', cmd)
