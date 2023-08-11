import re
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
    ):
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path, cram_path_string)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / cram_path_string),
        )

        assert sg.cram
        # ---- The job that we want to test
        return genotype(
            b=batch,
            sequencing_group_name=sg.id,
            cram_path=sg.cram,
            output_path=sg.gvcf,
            tmp_prefix=tmp_path,
        )

    def test_genotype_jobs_workflow_ref_fasta(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = 'workflow_overwritten_reference.fa'
        if isinstance(config.references['broad'], dict):
            config.references['broad']['ref_fasta'] = None
        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == 'HaplotypeCaller' or job.name == 'Postproc GVCF':
                assert re.search(r'workflow_overwritten_reference.fa', cmd)
                assert not re.search(r'hg38_reference.fa', cmd)

    def test_genotype_jobs_default_ref_fasta(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = None
        if isinstance(config.references['broad'], dict):
            config.references['broad']['ref_fasta'] = 'default_reference.fa'

        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == 'HaplotypeCaller' or job.name == 'Postproc GVCF':
                assert re.search(r'default_reference.fa', cmd)
                assert not re.search(r'hg38_reference.fa', cmd)

    def test_genotype_jobs_no_ref_fasta(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = None
        if isinstance(config.references['broad'], dict):
            config.references['broad']['ref_fasta'] = None

        with raises(ConfigError) as e:
            self._get_new_genotype_job(tmp_path, config)
        assert e.match(r'Failed to get reference fasta')

    def test_genotype_jobs_both_references_set(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = 'workflow_overwritten_reference.fa'
        if isinstance(config.references['broad'], dict):
            config.references['broad']['ref_fasta'] = 'default_reference.fa'

        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == 'HaplotypeCaller' or job.name == 'Postproc GVCF':
                assert re.search(r'workflow_overwritten_reference.fa', cmd)
                assert not re.search(r'default_references.fa', cmd)

    def test_genotype_jobs_with_default_scatter_count(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        # genotype_jobs = self._get_new_genotype_job(tmp_path, config)
        set_config(config, tmp_path / 'config.toml')

        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / 'test_genotype.cram'),
        )
        assert sg.cram
        # ---- The job that we want to test
        genotype_jobs = genotype(
            b=batch,
            sequencing_group_name=sg.id,
            cram_path=sg.cram,
            output_path=sg.gvcf,
            tmp_prefix=tmp_path,
        )

        # ---- Assertions
        expected_scatter_count = 50
        merge_job_name = f'Merge {expected_scatter_count} GVCFs'
        postproc_job_name = 'Postproc GVCF'
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for genome'
        assert (
            len(genotype_jobs) == expected_scatter_count + 3
        )  # +3 for making intervals, merging and postprocessing
        job_names = [job.name for job in genotype_jobs]
        assert make_intervals_job_name in job_names
        assert merge_job_name in job_names
        assert postproc_job_name in job_names

        haplotype_output_paths: list[str] = []
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == make_intervals_job_name:
                assert re.search(r'picard', cmd)
                assert re.search(r'IntervalListTools', cmd)
                assert re.search(r'SCATTER_COUNT=50', cmd)

                # Testing that outputs will be unique, sorted and internally divided
                assert re.search(r'SUBDIVISION_MODE=INTERVAL_SUBDIVISION', cmd)
                assert re.search(r'UNIQUE=true', cmd)
                assert re.search(r'SORT=true', cmd)

            elif job.name == merge_job_name:
                assert re.search(r'picard', cmd)
                assert re.search(r'MergeVcfs', cmd)
                # Validate the outputs from the haploytype caller jobs are the inputs to the merge job
                for haplotype_output_path in haplotype_output_paths:
                    assert re.search(re.escape(haplotype_output_path), cmd)

                # Validate the output is unique
                assert len(set(haplotype_output_paths)) == len(haplotype_output_paths)
                assert len(haplotype_output_paths) == expected_scatter_count

            elif job.name == postproc_job_name:
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

            else:
                # HaplotypeCaller jobs
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                assert re.search(r'--dragen-mode', cmd)
                pattern = r'-O\s+([^\\]+\.gz)'
                match = re.search(pattern, cmd)
                assert match
                haplotype_output_paths.append(match.group(1))

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

    def test_genotype_jobs_with_custom_scatter_count(self, tmp_path: Path):
        config = default_config()
        config.workflow.scatter_count_genotype = 10
        set_config(config, tmp_path / 'config.toml')

        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        expected_scatter_count = 10
        merge_job_name = f'Merge {expected_scatter_count} GVCFs'
        postproc_job_name = 'Postproc GVCF'
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for genome'
        assert (
            len(genotype_jobs) == expected_scatter_count + 3
        )  # +3 for making intervals, merging and postprocessing
        job_names = [job.name for job in genotype_jobs]
        assert make_intervals_job_name in job_names
        assert merge_job_name in job_names
        assert postproc_job_name in job_names

        haplotype_output_paths: list[str] = []
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == merge_job_name:
                assert re.search(r'picard', cmd)
                assert re.search(r'MergeVcfs', cmd)
                # Validate the outputs from the haploytype caller jobs are the inputs to the merge job
                for haplotype_output_path in haplotype_output_paths:
                    assert re.search(re.escape(haplotype_output_path), cmd)

                # Validate the output is unique
                assert len(set(haplotype_output_paths)) == len(haplotype_output_paths)
                assert len(haplotype_output_paths) == expected_scatter_count
            if job.name == 'HaplotypeCaller':
                # HaplotypeCaller jobs
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                pattern = r'-O\s+([^\\]+\.gz)'
                match = re.search(pattern, cmd)
                assert match
                haplotype_output_paths.append(match.group(1))

    def test_genotype_jobs_with_invalid_scatter_count(self, tmp_path: Path):
        config = default_config()
        config.workflow.scatter_count_genotype = 0
        set_config(config, tmp_path / 'config.toml')

        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        print(genotype_jobs)

        # ---- Assertions
        expected_scatter_count = 0
        merge_job_name = f'Merge {expected_scatter_count} GVCFs'
        postproc_job_name = 'Postproc GVCF'
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for genome'
        assert (
            len(genotype_jobs) == 2
        )  # We should only do the Postproc GVCF job and one haplotype caller job

        job_names = [job.name for job in genotype_jobs]
        assert make_intervals_job_name not in job_names
        assert merge_job_name not in job_names
        assert postproc_job_name in job_names

        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == postproc_job_name:
                assert re.search(r'gatk', cmd)
                assert re.search(r'ReblockGVCF', cmd)
            else:
                # HaplotypeCaller job
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                pattern = r'-O\s+([^\\]+\.gz)'
                match = re.search(pattern, cmd)
                assert match

    def test_genotype_default_reblock_gq_bands(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        postproc_job_name = 'Postproc GVCF'

        for job in genotype_jobs:
            if job.name == postproc_job_name:
                cmd = get_command_str(job)
                assert re.search(r'ReblockGVCF', cmd)
                assert re.search(
                    r'-floor-blocks -GQB 20 -GQB 30 -GQB 40',
                    cmd,
                )

    def test_genotype_custom_reblock_gq_bands(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.reblock_gq_bands = [10, 20, 50]
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        postproc_job_name = 'Postproc GVCF'
        for job in genotype_jobs:
            if job.name == postproc_job_name:
                cmd = get_command_str(job)
                assert re.search(r'ReblockGVCF', cmd)
                assert re.search(
                    r'-floor-blocks -GQB 10 -GQB 20 -GQB 50',
                    cmd,
                )

    def test_interval_tool_list_for_exome_with_defaults(self, tmp_path: Path):
        # ---- Test setup

        config = default_config()
        config.workflow.sequencing_type = 'exome'
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        expected_scatter_count = 50
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for exome'
        assert make_intervals_job_name in [job.name for job in genotype_jobs]
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == make_intervals_job_name:
                assert re.search(r'BREAK_BANDS_AT_MULTIPLES_OF=0 ', cmd)

    def test_interval_tool_list_for_genome_with_defaults(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.sequencing_type = 'genome'
        set_config(config, tmp_path / 'config.toml')
        genotype_jobs = self._get_new_genotype_job(tmp_path, config)

        # ---- Assertions
        expected_scatter_count = 50
        make_intervals_job_name = f'Make {expected_scatter_count} intervals for genome'
        assert make_intervals_job_name in [job.name for job in genotype_jobs]
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == make_intervals_job_name:
                assert re.search(r'BREAK_BANDS_AT_MULTIPLES_OF=100000 ', cmd)

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
