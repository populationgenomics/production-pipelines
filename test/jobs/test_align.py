"""
Testing Alignment Job
"""
from pathlib import Path

import pytest
from pytest_mock import MockFixture

from cpg_workflows.jobs.align import align
from cpg_workflows.resources import STANDARD


from test import set_config
from test.factories.config import create_config
from test.factories.batch import create_local_batch
from test.factories.sequencing_group import create_sequencing_group
from test.factories.alignment_input import (
    create_fastq_pair_input,
    create_fastq_pairs_input,
    create_cram_input,
    create_bam_input,
)


# ------------------------------------------------------------------------------------ #
# Tests
# ------------------------------------------------------------------------------------ #
class TestSingleFastqAlignment:
    def test_creates_one_alignment_job_and_one_markdupes_job(self, tmp_path: Path):
        config = create_config(sequencing_type="genome")
        set_config(config, tmp_path / 'config.toml')

        batch = create_local_batch(tmp_path)
        sequencing_group = create_sequencing_group(
            alignment_input=create_fastq_pair_input(location=tmp_path),
        )

        jobs = align(b=batch, sequencing_group=sequencing_group)

        assert len(jobs) == 2
        assert "Align" in str(jobs[0].name)
        assert "MarkDuplicates" in str(jobs[1].name)


# class TestMultipleFastqlignment:
#     def test_creates_one_alignment_job_for_each_pair(self, tmp_path: Path):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')
#         batch = create_local_batch(tmp_path)

#         sequencing_group = create_sequencing_group(
#             alignment_input_type='fastq-pairs',
#             number_of_fastq_pairs=4,
#         )

#         jobs = align(b=batch, sequencing_group=sequencing_group)

#         align_jobs = [job for job in jobs if "Align" in str(job.name)]
#         assert len(align_jobs) == 4

#     def test_creates_a_merge_job(self, tmp_path: Path):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')
#         batch = create_local_batch(tmp_path)

#         sequencing_group = create_sequencing_group(
#             alignment_input_type='fastq-pairs', number_of_fastq_pairs=4
#         )

#         jobs = align(b=batch, sequencing_group=sequencing_group)

#         merge_job = [job for job in jobs if "Merge" in str(job.name)]
#         assert len(merge_job) == 1

#     def test_creates_a_markdupes_job(self, tmp_path: Path):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')
#         batch = create_local_batch(tmp_path)

#         sequencing_group = create_sequencing_group(
#             alignment_input_type='fastq-pairs', number_of_fastq_pairs=4
#         )

#         jobs = align(b=batch, sequencing_group=sequencing_group)

#         markdupes_job = [job for job in jobs if "MarkDuplicates" in str(job.name)]
#         assert len(markdupes_job) == 1


# class TestRealignmentFromCram:
#     pass


# class TestRealignmentFromBam:
#     pass


# class TestJobMetadata:
#     def test_adds_extra_label_to_align_job_name(self, tmp_path: Path):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')

#         sequencing_group = create_sequencing_group()
#         batch = create_local_batch(tmp_path)

#         jobs = align(b=batch, sequencing_group=sequencing_group, extra_label="MyTest")
#         align_jobs = [job for job in jobs if "Align" in str(job.name)]
#         for job in align_jobs:
#             assert 'MyTest' in str(job.name)


# class TestResourceAquisition:
#     def test_align_job_requests_resources_from_standard_machine_type(
#         self, mocker: MockFixture, tmp_path: Path
#     ):
#         config = create_config(workflow_resources={"storage_gb": 100})
#         set_config(config, tmp_path / 'config.toml')

#         sequencing_group = create_sequencing_group()
#         batch = create_local_batch(tmp_path)

#         spy = mocker.spy(STANDARD, "set_resources")
#         jobs = align(b=batch, sequencing_group=sequencing_group, requested_nthreads=4)
#         spy.assert_called_with(jobs[0], nthreads=4, storage_gb=100)

#     # def test_using_biobambam_saves_markdup_metrics():
#     #     # ... bunch of code up that sets spy write_output

#     #     import re

#     #     re.match(r"M=.*/markdup_metrics")
#     #     spy.assert_called_with(jobs[1].markdup_metrics, out_markdup_metrics_path)

#     def test_align_job_allocates_additional_memory_for_genome_alignments(
#         self, tmp_path: Path
#     ):
#         # Get genome storage
#         config = create_config(sequencing_type="genome")
#         set_config(config, tmp_path / 'genome_config.toml')

#         batch = create_local_batch(tmp_path)
#         sequencing_group = create_sequencing_group(sequencing_type="genome")

#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         genome_storage = float(str(jobs[0]._storage).replace('G', ''))

#         # Get exome storage
#         config = create_config(sequencing_type="exome")
#         sequencing_group = create_sequencing_group(sequencing_type="exome")
#         set_config(config, tmp_path / 'exome_config.toml')
#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         exome_storage = float(str(jobs[0]._storage).replace('G', ''))

#         assert genome_storage > exome_storage

#     @pytest.mark.parametrize(
#         "input_type",
#         [
#             "bam-with-index",
#             "fastq-pair",
#             "fastq-pairs",
#         ],
#     )
#     def test_align_job_allocates_less_storage_for_indexed_crams(
#         self, tmp_path: Path, input_type: AlignmentInputType
#     ):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')

#         batch = create_local_batch(tmp_path)

#         # Get CRAM storage
#         sequencing_group = create_sequencing_group(
#             alignment_input_type="cram-with-index"
#         )
#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         cram_storage = float(str(jobs[0]._storage).replace('G', ''))

#         # Get other storage
#         sequencing_group = create_sequencing_group(alignment_input_type=input_type)
#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         storage = float(str(jobs[0]._storage).replace('G', ''))
#         assert (
#             storage > cram_storage
#         ), f"'{input_type}' storage should be greater than cram storage"

#     @pytest.mark.parametrize(
#         "with_index,without_index",
#         [
#             ("bam-with-index", "bam-without-index"),
#             ("cram-with-index", "cram-without-index"),
#         ],
#     )
#     def test_align_job_allocates_more_storage_for_unindexed_crams_bams(
#         self,
#         tmp_path: Path,
#         with_index: AlignmentInputType,
#         without_index: AlignmentInputType,
#     ):
#         config = create_config()
#         set_config(config, tmp_path / 'config.toml')

#         batch = create_local_batch(tmp_path)
#         sequencing_group = create_sequencing_group(alignment_input_type=with_index)
#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         storage_with_index = float(str(jobs[0]._storage).replace('G', ''))

#         sequencing_group = create_sequencing_group(alignment_input_type=without_index)
#         jobs = align(b=batch, sequencing_group=sequencing_group)
#         storage_without_index = float(str(jobs[0]._storage).replace('G', ''))

#         assert (
#             storage_without_index > storage_with_index
#         ), f"'{without_index}' storage should be greater than '{with_index}' storage"
