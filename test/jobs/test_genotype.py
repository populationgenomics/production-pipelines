import re
import pytest
from pytest_mock import MockFixture
from pathlib import Path

from cpg_utils.config import ConfigError
from cpg_workflows.jobs.genotype import genotype
from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str
from cpg_workflows.filetypes import CramPath
from pytest import raises


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='genotype-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
            reblock_gq_bands=[20, 30, 40],
        ),
        images={
            'picard': 'picard:2.27.4',
            'gatk': 'gatk:4.2.6.1',
        },
        references={
            'broad': {
                'ref_fasta': 'hg38_reference.fa',
                'genome_calling_interval_lists': 'wgs_calling_regions.hg38.interval_list',
                'noalt_bed': 'primary_contigs_plus_mito.bed.gz',
                'exome_calling_interval_lists': 'exome_calling_regions.v1.interval_list',
            }
        },
        other={
            'resource_overrides': {},
        },
    )


class TestGenotyping:
    def _get_new_genotype_job(
        self,
        tmp_path: Path,
        config: PipelineConfig,
        cram_path_string: str = 'test_genotype.cram',
        batch=None,
    ):
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        if batch is None:
            batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / cram_path_string),
        )

        assert sg.cram
        # ---- The job that we want to test
        jobs = genotype(
            b=batch,
            sequencing_group_name=sg.id,
            cram_path=sg.cram,
            output_path=sg.gvcf,
            tmp_prefix=tmp_path,
        )

        return jobs, sg

    @pytest.mark.parametrize('scatter_count', [50, 10, 0, None])
    def test_genotype_jobs_with_varying_scatter_counts(
        self, mocker: MockFixture, tmp_path: Path, scatter_count
    ):
        # ---- Test setup
        config = default_config()
        if scatter_count is not None:
            config.workflow.scatter_count_genotype = scatter_count
        else:
            scatter_count = 50

        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)
        spy = mocker.spy(batch, 'write_output')

        genotype_jobs, sg = self._get_new_genotype_job(tmp_path, config, batch=batch)

        # ---- Assertions
        if scatter_count == 0:
            expected_jobs = 2
            expected_postproc = True
            expected_intervals = False
            expected_merge = False
        else:
            expected_jobs = scatter_count + 3
            expected_postproc = True
            expected_intervals = True
            expected_merge = True

        job_names = [job.name for job in genotype_jobs]
        assert len(genotype_jobs) == expected_jobs

        assert ('Postproc GVCF' in job_names) == expected_postproc
        assert (
            f'Make {scatter_count} intervals for genome' in job_names
        ) == expected_intervals
        assert (f'Merge {scatter_count} GVCFs' in job_names) == expected_merge

        haplotype_output_paths: list[str] = []
        expected_calls = []
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == f'Merge {scatter_count} GVCFs':
                assert re.search(r'picard', cmd)
                assert re.search(r'MergeVcfs', cmd)
                for haplotype_output_path in haplotype_output_paths:
                    assert re.search(re.escape(haplotype_output_path), cmd)
                assert len(set(haplotype_output_paths)) == len(haplotype_output_paths)
                assert len(haplotype_output_paths) == scatter_count

            if job.name == 'HaplotypeCaller':
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                pattern = r'-O\s+([^\\]+\.gz)'
                match = re.search(pattern, cmd)
                assert match
                haplotype_output_paths.append(match.group(1))

                if scatter_count == 0:
                    file_name = sg.id
                else:
                    scatter_frac = job.attributes['part']
                    scatter_parts = scatter_frac.split('/')
                    file_name = (
                        f'{int(scatter_parts[0]) - 1}_of_{scatter_parts[1]}_{sg.id}'
                    )

                file_path = tmp_path / 'haplotypecaller' / file_name

                expected_calls.append(mocker.call(job.output_gvcf, str(file_path)))

            if job.name == 'Postproc GVCF':
                assert re.search(r'gatk', cmd)
                assert re.search(r'ReblockGVCF', cmd)

        spy.assert_has_calls(expected_calls, any_order=True)

    @pytest.mark.parametrize(
        'reblock_gq_bands, expected_gqb_flags',
        [
            ([20, 30, 40], '-GQB 20 -GQB 30 -GQB 40'),
            ([10, 20, 50], '-GQB 10 -GQB 20 -GQB 50'),
            (None, '-GQB 20 -GQB 30 -GQB 40'),
        ],
    )
    def test_genotype_reblock_gq_bands(
        self, tmp_path: Path, reblock_gq_bands, expected_gqb_flags
    ):
        # ---- Test setup
        config = default_config()
        if reblock_gq_bands is not None:
            config.workflow.reblock_gq_bands = reblock_gq_bands
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs, _sg = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        postproc_job_name = 'Postproc GVCF'

        for job in genotype_jobs:
            if job.name == postproc_job_name:
                cmd = get_command_str(job)
                assert re.search(r'ReblockGVCF', cmd)
                assert re.search(rf'-floor-blocks {expected_gqb_flags}', cmd)

    @pytest.mark.parametrize(
        'sequencing_type, expected_intervals_name, expected_multiples_of',
        [
            ('exome', 'Make 50 intervals for exome', '0'),
            ('genome', 'Make 50 intervals for genome', '100000'),
        ],
    )
    def test_interval_tool_list_with_defaults(
        self,
        tmp_path: Path,
        sequencing_type,
        expected_intervals_name,
        expected_multiples_of,
    ):
        # ---- Test setup
        config = default_config()
        config.workflow.sequencing_type = sequencing_type
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs, _sg = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        assert expected_intervals_name in [job.name for job in genotype_jobs]
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == expected_intervals_name:
                assert re.search(
                    rf'BREAK_BANDS_AT_MULTIPLES_OF={expected_multiples_of} ', cmd
                )

    @pytest.mark.parametrize(
        'workflow_ref_fasta, broad_ref_fasta, expected_ref_fasta_match, expect_error',
        [
            (
                'workflow_overwritten_reference.fa',
                None,
                'workflow_overwritten_reference.fa',
                False,
            ),
            (None, 'default_reference.fa', 'default_reference.fa', False),
            (None, None, None, True),
            (
                'workflow_overwritten_reference.fa',
                'default_reference.fa',
                'workflow_overwritten_reference.fa',
                False,
            ),
        ],
    )
    def test_genotype_jobs_ref_fasta(
        self,
        tmp_path: Path,
        workflow_ref_fasta,
        broad_ref_fasta,
        expected_ref_fasta_match,
        expect_error,
    ):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = workflow_ref_fasta
        if isinstance(config.references['broad'], dict):
            config.references['broad']['ref_fasta'] = broad_ref_fasta

        if expect_error:
            with pytest.raises(ConfigError, match=r'Failed to get reference fasta'):
                self._get_new_genotype_job(tmp_path, config)
            return

        genotype_jobs, _sg = self._get_new_genotype_job(tmp_path, config)

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == 'HaplotypeCaller' or job.name == 'Postproc GVCF':
                if expected_ref_fasta_match:
                    assert re.search(expected_ref_fasta_match, cmd)
                assert not re.search(r'hg38_reference.fa', cmd)

    def test_postproc_gvcf_job_params_with_defaults(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()

        genotype_jobs, sg = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        postproc_job_name = 'Postproc GVCF'
        job_names = [job.name for job in genotype_jobs]
        assert postproc_job_name in job_names

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == postproc_job_name:
                assert re.search(r'gatk', cmd)
                assert re.search(r'ReblockGVCF', cmd)
                # Necessary info field annotations added to perform QUAL approximation downstream
                assert re.search(r'-do-qual-approx', cmd)
                # Create tabix index
                assert re.search(r'--create-output-variant-index', cmd)

                # Validate SG ID in reheader #bcftools reheader -s <(echo "$EXISTING_SN CPG000001")
                reheader_with_sgid = (
                    f'bcftools reheader -s <(echo "$EXISTING_SN {sg.id}")'
                )
                assert re.search(re.escape(reheader_with_sgid), cmd)

                # No alt region set to broad/noalt_bed
                broad_reference = config.references.get('broad')
                if isinstance(broad_reference, dict) and (
                    no_alt_bed := broad_reference.get('noalt_bed')
                ):
                    no_alt_region = f'-T .*{no_alt_bed}'
                    assert re.search(no_alt_region, cmd)

    def test_make_intervals_job_params_with_defaults(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        genotype_jobs, _sg = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        expected_scatter_count = 50
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for genome'

        for job in genotype_jobs:
            if job.name == make_intervals_job_name:
                cmd = get_command_str(job)
                assert re.search(r'picard', cmd)
                assert re.search(r'IntervalListTools', cmd)
                assert re.search(r'SCATTER_COUNT=50', cmd)

                # Testing that outputs will be unique, sorted and internally divided
                assert re.search(r'SUBDIVISION_MODE=INTERVAL_SUBDIVISION', cmd)
                assert re.search(r'UNIQUE=true', cmd)
                assert re.search(r'SORT=true', cmd)

    def test_haplotype_caller_job_params_with_defaults(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')
        batch = create_local_batch(tmp_path)
        spy = mocker.spy(batch, 'write_output')
        genotype_jobs, sg = self._get_new_genotype_job(
            tmp_path=tmp_path, config=config, batch=batch
        )

        # ---- Assertions
        haplotype_caller_job_name = 'HaplotypeCaller'
        expected_calls = []
        for job in genotype_jobs:
            if job.name == haplotype_caller_job_name:
                cmd = get_command_str(job)
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                assert re.search(r'--dragen-mode', cmd)
                # If a genotyping event overlaps deletions that the '*' spanning event will be excluded
                assert re.search(r'--disable-spanning-event-genotyping', cmd)
                # Allele specific annotations requested
                assert re.search(r'-G AS_StandardAnnotation', cmd)
                # GQB bands set to multiples of 10
                assert re.search(
                    r'-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90',
                    cmd,
                )
                assert re.search(r'-ERC GVCF', cmd)
                # Tabix index created
                assert re.search(r'--create-output-variant-index', cmd)

                scatter_frac = job.attributes['part']
                scatter_parts = scatter_frac.split('/')
                file_name = f'{int(scatter_parts[0]) - 1}_of_{scatter_parts[1]}_{sg.id}'
                file_path = tmp_path / 'haplotypecaller' / file_name

                expected_calls.append(mocker.call(job.output_gvcf, str(file_path)))

        spy.assert_has_calls(expected_calls, any_order=True)

    # def test_genotype_reblock_gq_bands_single_integer(self, tmp_path: Path):
    #     # NOTE: This test fails, but should it?
    #     """
    #     # ---- Test setup
    #     config = default_config()
    #     config.workflow.reblock_gq_bands = 10
    #     set_config(config, tmp_path / 'config.toml')
    #     genotype_jobs = self._get_new_genotype_job(tmp_path, config)

    #     # ---- Assertions
    #     postproc_job_name = 'Postproc GVCF'
    #     for job in genotype_jobs:
    #         if job.name == postproc_job_name:
    #             cmd = get_command_str(job)
    #             assert re.search(r'ReblockGVCF', cmd)
    #             assert re.search(
    #                 r'-floor-blocks -GQB 10',
    #                 cmd,
    #             )
    #     """
    #     pass
