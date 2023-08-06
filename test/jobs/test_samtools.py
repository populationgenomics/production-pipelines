import os
import re

from cpg_utils import Path
from cpg_utils.hail_batch import image_path
from pytest_mock import MockFixture

from cpg_workflows import utils
from cpg_workflows.batch import Batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs.samtools import samtools_stats

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='samtools-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'samtools': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                }
            },
        },
    )


def setup_test(tmp_path: Path):
    config = default_config()
    set_config(config, tmp_path / 'config.toml')

    cram_pth = create_cram_input(
        location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
    )

    batch = create_local_batch(tmp_path)

    return config, cram_pth, batch


class TestSamtoolsRun:
    def test_samtools_stats_creates_a_job(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )

        # ---- Assertions
        assert (
            j is not None
        ), 'The samtools_stats function did not create a job. Check if overwrite=False in samtools_stats call'

    def test_creates_one_job(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # ---- Assertions
        assert (
            len(batch.select_jobs('samtools stats')) == 1
        ), "Unexpected number of 'samtools stats' jobs in batch list, should be just 1 job"

    def test_can_reuse_existing_output_path(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j1 = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        j2 = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=False,
        )

        # ---- Assertions
        assert (
            j2 is None
        ), 'New directory has been created for second job when it should reuse directory \
        of j1. Try setting overwrite to False in subsequent jobs.'

    def test_will_return_None_if_path_already_exists(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_stats_file'
        output.touch()

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=False,
        )

        # --- Assertions
        assert j is None

    def test_will_create_job_if_path_already_exists_and_overwrite_true(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)
        output = tmp_path / 'output_stats_file'
        output.touch()

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=output,
            overwrite=True,
        )

        # --- Assertions
        assert j is not None

    def test_pass_file_that_doesnt_exist(self, mocker: MockFixture, tmp_path: Path):
        # can_reuse() executes all(exists()) which checks whether all files in the path exist
        # if any of the files in the paths do not exist can_reuse() returns False and a job is created

        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=(tmp_path / 'output_stats_file'),
            job_attrs=None,
            overwrite=False,
        )

        # ---- Assertions
        assert j is not None

    def test_can_overwrite_existing_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        # Create a spy on cpg_workflows.utils.can_reuse
        # TODO: why is spy not working?
        # spy_can_reuse = mocker.spy(utils, 'can_reuse')

        initial_cram_pth = create_cram_input(
            location=tmp_path,
            prefix='test1',
            index=True,
            reference_assembly='GRCh38.fa',
            create=True,
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j1 = samtools_stats(
            b=batch,
            cram_path=initial_cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # ---- Assert initial .cram & .crai paths are created and point to correct directory
        assert initial_cram_pth.path == (tmp_path / 'test1.cram')
        assert initial_cram_pth.index_path == (tmp_path / 'test1.cram.crai')

        second_cram_pth = create_cram_input(
            location=tmp_path,
            prefix='test1',
            index=True,
            reference_assembly='GRCh38.fa',
            create=True,
        )

        j2 = samtools_stats(
            b=batch,
            cram_path=second_cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # ---- Assertions
        # assert spy_can_reuse.assert_called_once_with(tmp_path, True)
        assert second_cram_pth.path == initial_cram_pth.path
        assert second_cram_pth.index_path == initial_cram_pth.index_path

    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j_default_attrs = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        j_supplied_attrs = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs={'test_tool': 'test_samtools'},
            overwrite=True,
        )

        assert j_default_attrs.attributes == {'tool': 'samtools'}
        assert j_supplied_attrs.attributes == {
            'test_tool': 'test_samtools',
            'tool': 'samtools',
        }

    def test_uses_samtools_image_specified_in_config(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        assert j._image == config.images['samtools']

    # TODO: reference is '__RESOURCE_FILE__0' how do i check?
    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        # Not sure how to check what reference file is being used. \
        # The resource group inside the job 'j' is not obvious in its naming

        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        pass

    # TODO: same as above
    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section():
        pass

    # TODO: what is meant by this test?
    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        print()

    # Can't test this properly as no output is being generated
    def test_samtools_writes_to_resource_file_named_output_stats(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # ---- Assertions
        cmd = get_command_str(j)
        match = re.search(r'\$CRAM > (\S+)', cmd)
        # Retrieve the j.output_stats from the regex match
        actual_output_stats = match.group(1)
        assert (
            actual_output_stats == j.output_stats
        ), f'Expected: {j.output_stats}, Got: {actual_output_stats}'

    # TODO: What's the difference between this test and testing whether output exists?
    def test_batch_writes_samtools_stats_file_to_output_path():
        pass

    # TODO: Check if implemented spy object properly
    def test_check_batch_command_executed(self, mocker: MockFixture, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # Set up spy object to track running of batch jobs
        spy_batch_run = mocker.spy(Batch, 'run')

        # ---- Assertions
        spy_batch_run.assert_called_once

    def test_check_output_files_exist(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        set_config(config, tmp_path / 'config.toml')

        cram_pth = create_cram_input(
            location=tmp_path, prefix='test', index=True, reference_assembly='GRCh38.fa'
        )

        batch = create_local_batch(tmp_path)

        # ---- The jobs we want to test
        j = samtools_stats(
            b=batch,
            cram_path=cram_pth,
            out_samtools_stats_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )
        # ---- Assertions
        assert (
            tmp_path / f'{j.output_stats}'
        ).exists(), f'No output file named {j.output_stats}'
