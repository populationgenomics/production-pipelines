import re
from pathlib import Path

from typing import Literal
from cpg_utils.config import ConfigError
from cpg_workflows.jobs.happy import happy
from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str
from cpg_workflows.filetypes import CramPath, GvcfPath
from pytest import raises
import pytest
from pytest_mock import MockFixture


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='happy-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'hap-py': 'hap-py:0.3.15',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'genome_calling_interval_lists': 'wgs_calling_regions.hg38.interval_list',
                'noalt_bed': 'primary_contigs_plus_mito.bed.gz',
                'exome_calling_interval_lists': 'exome_calling_regions.v1.interval_list',
                'genome_evaluation_interval_lists': 'genome_evaluation_regions.v1.interval_list',
                'exome_evaluation_interval_lists': 'exome_evaluation_regions.v1.interval_list',
            },
            'sample1': {
                'truth_vcf': 'sample_one_truth.vcf.gz',
                'regions_bed': 'sample_one_regions.bed',
            },
            'sample2': {
                'truth_vcf': 'sample_two_truth.vcf.gz',
                'regions_bed': 'sample_two_regions.bed',
            },
        },
        other={
            'resource_overrides': {},
            'validation': {'sample_map': {'SAMPLE1': 'sample1'}},
        },
    )


class TestHappy:
    def _get_new_happy_job(
        self,
        tmp_path: Path,
        config: PipelineConfig,
        is_gvcf: bool = True,
    ):
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path, 'test-batch')

        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / 'sample_one.cram'),
            # Can't use tmp_path for gvcf, because it expects a gs:// path
            gvcf=GvcfPath('gs://test-project/gvcf/sample_one_genotype.g.vcf.gz'),
        )

        assert sg.gvcf

        if is_gvcf:
            happy_job = happy(
                sequencing_group=sg,
                b=batch,
                vcf_or_gvcf=GvcfPath(sg.gvcf).resource_group(batch),
                is_gvcf=is_gvcf,
                job_attrs={},
            )
        else:
            vcf_path = tmp_path / 'test_vcf.vcf.gz'
            happy_job = happy(
                sequencing_group=sg,
                b=batch,
                vcf_or_gvcf=batch.read_input_group(
                    **{
                        'vcf.gz': str(vcf_path),
                        'vcf.gz.tbi': str(vcf_path) + '.tbi',
                    }
                ),
                is_gvcf=is_gvcf,
                job_attrs={},
            )

        return happy_job

    def test_empty_validation_config_returns_no_job(self, tmp_path: Path):
        config = default_config()
        config.other['validation'] = {}
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job is None

    def test_sample_missing_from_validation_config_returns_no_job(self, tmp_path: Path):
        config = default_config()
        config.other['validation'] = {'sample_map': {'SAMPLE2': 'sample2'}}
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job is None

    def test_genome_sequencing_type_gvcf(self, tmp_path: Path):
        config = default_config()
        config.workflow.sequencing_type = 'genome'
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job
        assert happy_job.name == 'hap.py (GVCF)'
        cmd = get_command_str(happy_job)
        assert re.search(r'--restrict-regions', cmd)
        assert not re.search(r'--target-regions', cmd)

    def test_genome_sequencing_type_vcf(self, tmp_path: Path):
        config = default_config()
        config.workflow.sequencing_type = 'genome'
        happy_job = self._get_new_happy_job(tmp_path, config, is_gvcf=False)
        assert happy_job
        assert happy_job.name == 'hap.py (VCF)'
        cmd = get_command_str(happy_job)
        assert re.search(r'--restrict-regions', cmd)
        assert not re.search(r'--target-regions', cmd)

    def test_exome_sequencing_type_gvcf(self, tmp_path: Path):
        config = default_config()
        config.workflow.sequencing_type = 'exome'
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job
        assert happy_job.name == 'hap.py (GVCF)'
        cmd = get_command_str(happy_job)
        assert re.search(r'--target-regions', cmd)
        assert not re.search(r'--restrict-regions', cmd)

    def test_exome_sequencing_type_vcf(self, tmp_path: Path):
        config = default_config()
        config.workflow.sequencing_type = 'exome'
        happy_job = self._get_new_happy_job(tmp_path, config, is_gvcf=False)
        assert happy_job
        assert happy_job.name == 'hap.py (VCF)'
        cmd = get_command_str(happy_job)
        assert re.search(r'--target-regions', cmd)
        assert not re.search(r'--restrict-regions', cmd)

    def test_valid_sample_but_excluded_from_config_returns_no_job(self, tmp_path: Path):
        config = default_config()
        config.other['validation'] = {'sample_map': {'SAMPLE2': 'sample2'}}
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job is None

    def test_correct_sample_is_selected_from_validation_config(self, tmp_path: Path):
        config = default_config()
        config.other['validation'] = {
            'sample_map': {'SAMPLE2': 'sample2', 'SAMPLE1': 'sample1'}
        }
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job.name == 'hap.py (GVCF)'
        cmd = get_command_str(happy_job)
        assert re.search(r'sample_one_genotype.g.vcf.gz', cmd)
        assert re.search(r'sample_one_truth.vcf.gz', cmd)
        assert re.search(r'sample_one_regions.bed', cmd)
        assert not re.search(r'sample_two_truth.vcf', cmd)
        assert not re.search(r'sample_two_regions.bed', cmd)

    def test_happy_params_with_gvcf_input(self, tmp_path: Path):
        config = default_config()
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job.name == 'hap.py (GVCF)'
        assert happy_job._image == config.images['hap-py']
        cmd = get_command_str(happy_job)
        assert re.search(r'hap.py', cmd)
        assert re.search(r'--convert-gvcf-to-vcf', cmd)
        assert re.search(r'--filter-nonref', cmd)
        assert not re.search(r'bcftools view', cmd)
        assert not re.search(r'bcftools annotate', cmd)
        assert not re.search(r'bcftools index', cmd)
        assert re.search(r'sample_one_genotype.g.vcf.gz', cmd)

    def test_happy_params_with_vcf_input(self, tmp_path: Path):
        config = default_config()
        happy_job = self._get_new_happy_job(tmp_path, config, is_gvcf=False)
        assert happy_job.name == 'hap.py (VCF)'
        assert happy_job._image == config.images['hap-py']
        cmd = get_command_str(happy_job)
        assert re.search(r'hap.py', cmd)
        assert not re.search(r'--convert-gvcf-to-vcf', cmd)
        assert not re.search(r'--filter-nonref', cmd)
        assert re.search(r'bcftools view', cmd)
        assert re.search(r'bcftools annotate', cmd)
        assert re.search(r'bcftools index', cmd)
        assert re.search(r'test_vcf.vcf.gz', cmd)

    def test_skips_happy_if_validation_not_set(self, tmp_path: Path):
        config = default_config()
        config.other['validation'] = {}
        happy_job = self._get_new_happy_job(tmp_path, config)
        assert happy_job is None

    def test_invalid_references_raises_config_error(self, tmp_path: Path):
        config = default_config()
        config.references['broad'] = {}
        with raises(ConfigError):
            self._get_new_happy_job(tmp_path, config)

    @pytest.mark.parametrize('sequencing_type', ['genome', 'exome'])
    def test_write_output_is_correctly_called_when_output_path_specified(
        self,
        mocker: MockFixture,
        tmp_path: Path,
        sequencing_type: Literal['genome', 'exome'],
    ):
        config = default_config()
        config.workflow.sequencing_type = sequencing_type
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path, 'test_batch')

        spy = mocker.spy(batch, 'write_output')
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / 'sample_one.cram'),
            gvcf=GvcfPath('gs://test-project/gvcf/sample_one_genotype.g.vcf.gz'),
        )
        output_path = tmp_path / 'test_output.csv'
        happy_job = happy(
            sequencing_group=sg,
            b=batch,
            vcf_or_gvcf=GvcfPath(sg.gvcf).resource_group(batch),
            is_gvcf=True,
            job_attrs={},
            output_path=output_path,
        )

        spy.assert_called_once_with(happy_job.summary_csv, str(output_path))

    @pytest.mark.parametrize('sequencing_type', ['genome', 'exome'])
    @pytest.mark.parametrize(
        'external_id, sample_map',
        [('SAMPLE1', 'sample_one'), ('SAMPLE2', 'sample_two')],
    )
    def test_read_input_is_correctly_called(
        self,
        mocker: MockFixture,
        tmp_path: Path,
        sequencing_type: Literal['genome', 'exome'],
        external_id: str,
        sample_map: str,
    ):
        config = default_config()
        config.workflow.sequencing_type = sequencing_type
        config.other['validation']['sample_map'] = {external_id: external_id.lower()}
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path, 'test_batch')

        spy = mocker.spy(batch, 'read_input')

        sg = create_sequencing_group(
            external_id=external_id,
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / 'sample_one.cram'),
            gvcf=GvcfPath('gs://test-project/gvcf/sample_one_genotype.g.vcf.gz'),
        )
        happy_job = happy(
            sequencing_group=sg,
            b=batch,
            vcf_or_gvcf=GvcfPath(sg.gvcf).resource_group(batch),
            is_gvcf=True,
            job_attrs={},
        )

        assert happy_job

        spy.assert_has_calls(
            [
                mocker.call(f'{sequencing_type}_evaluation_regions.v1.interval_list'),
                mocker.call(f'{sample_map}_truth.vcf.gz'),
                mocker.call(f'{sample_map}_regions.bed'),
            ]
        )
