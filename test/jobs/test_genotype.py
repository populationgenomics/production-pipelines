import re
from pathlib import Path

from cpg_workflows.jobs.genotype import genotype

from .. import set_config
from ..factories.alignment_input import create_fastq_pairs_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str
from cpg_workflows.filetypes import CramPath


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='genotype-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
            reblock_gq_bands=[20, 30, 50],
        ),
        images={
            'picard': 'picard:2.27.4',
            'gatk': 'gatk:4.2.6.1',
        },
        other={
            'resource_overrides': {},
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'genome_calling_interval_lists': 'wgs_calling_regions.hg38.interval_list',
                    'noalt_bed': 'primary_contigs_plus_mito.bed.gz',
                }
            },
        },
    )


class TestGenotyping:
    def test_genotype_jobs_with_default_scatter_count(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
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

        # TODO: Add a list of params we could check
        haplotype_output_paths: list[str] = []
        for job in genotype_jobs:
            cmd = get_command_str(job)
            if job.name == make_intervals_job_name:
                # TODO: Validate input and output
                print(cmd)
                assert re.search(r'picard', cmd)
                assert re.search(r'IntervalListTools', cmd)
                assert re.search(r'SCATTER_COUNT=50', cmd)

            elif job.name == merge_job_name:
                assert re.search(r'picard', cmd)
                assert re.search(r'MergeVcfs', cmd)
                # Validate the outputs from the haploytype caller jobs are the inputs to the merge job
                for haplotype_output_path in haplotype_output_paths:
                    assert re.search(re.escape(haplotype_output_path), cmd)

            elif job.name == postproc_job_name:
                # TODO Validate input and output
                print(cmd)
                assert re.search(r'gatk', cmd)
                assert re.search(r'ReblockGVCF', cmd)
            else:
                # HaplotypeCaller jobs
                # TODO: Validate input and output
                assert re.search(r'gatk', cmd)
                assert re.search(r'HaplotypeCaller', cmd)
                pattern = r'-O\s+([^\\]+\.gz)'
                match = re.search(pattern, cmd)
                assert match
                haplotype_output_paths.append(match.group(1))
                print(cmd)

    def test_genotype_jobs_with_custom_scatter_count(self, tmp_path: Path):
        config = default_config()
        config.workflow.scatter_count_genotype = 10
        set_config(config, tmp_path / 'config.toml')
        dataset_id = config.workflow.dataset
        batch = create_local_batch(tmp_path)
        sg = create_sequencing_group(
            dataset=dataset_id,
            sequencing_type=config.workflow.sequencing_type,
            alignment_input=create_fastq_pairs_input(location=tmp_path, n=1),
            cram=CramPath(tmp_path / 'test_genotype2.cram'),
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
