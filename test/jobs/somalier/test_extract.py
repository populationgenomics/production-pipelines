import re

import pytest
from pytest_mock import MockFixture

from cpg_utils import Path
from cpg_workflows.jobs.somalier import extract

from ... import set_config
from ...factories.alignment_input import create_cram_input
from ...factories.batch import create_local_batch
from ...factories.config import PipelineConfig, WorkflowConfig
from ..helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='somalier-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'somalier': 'test_image',
            'cpg_workflows': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                },
                'somalier_sites': 'test_somalier_sites',
            },
            'cramqc': {'num_pcs': '4'},
        },
    )


def setup_test(
    tmp_path: Path,
    config: PipelineConfig | None = None,
    index: bool = True,
):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    cram_pth = create_cram_input(location=tmp_path, index=index)
    batch = create_local_batch(tmp_path)

    return config, cram_pth, batch


class TestSomalierExtract:
    def test_creates_one_job(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        _ = extract(b=batch, cram_path=cram_pth, out_somalier_path=tmp_path)

        assert (
            len(batch.select_jobs('Somalier extract')) == 1
        ), "Unexpected number of 'Somalier' jobs in batch list, should be just 1 job"

    def test_will_return_none_if_output_file_exists_and_overwrite_is_false(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=False,
        )

        assert j is None, 'Job was created when output file should have been reused'

    def test_will_create_one_job_if_output_file_exists_and_overwrite_is_true(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=True,
        )

        assert len(batch.select_jobs('Somalier extract')) == 1
        assert j is not None

    @pytest.mark.parametrize(
        'job_attrs, expected_attrs',
        [
            (None, {'tool': 'somalier'}),
            (
                {'test_tool': 'test_somalier'},
                {'test_tool': 'test_somalier', 'tool': 'somalier'},
            ),
        ],
    )
    def test_sets_supplied_job_attrs_or_sets_default_attrs_if_attrs_not_supplied(
        self,
        tmp_path: Path,
        job_attrs,
        expected_attrs,
    ):
        _, cram_pth, batch = setup_test(tmp_path)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
            job_attrs=job_attrs,
        )

        assert j is not None
        assert j.attributes == expected_attrs

    def test_uses_image_specified_in_config(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        assert j is not None
        assert j._image == config.images['somalier']

    def test_sets_sites_location_with_name_of_sites_file_in_config(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        cmd = get_command_str(j)
        sites = config.other['references']['somalier_sites']
        assert j is not None
        assert re.search(fr'SITES=\$BATCH_TMPDIR/sites/{sites}', cmd)
        assert re.search(fr'--sites \${{BATCH_TMPDIR}}/inputs/\w+/{sites}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(self, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        cmd = get_command_str(j)
        # This test ensures the script uses a fail-safe copy operation when working with CRAM files and their indices.
        # Such retries are common for network calls like copying files from GCS to handle unexpected issues.
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        config = default_config()
        config.workflow.ref_fasta = 'test_workflow_ref.fa'
        _, cram_pth, batch = setup_test(tmp_path, config)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(self, tmp_path: Path):
        config, cram_pth, batch = setup_test(tmp_path)

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        cmd = get_command_str(j)
        ref_file = config.other['references']['broad']['ref_fasta']
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_batch_writes_output_file_to_output_path(self, mocker: MockFixture, tmp_path: Path):
        _, cram_pth, batch = setup_test(tmp_path)

        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=out_path,
        )

        assert j is not None
        spy.assert_called_with(j.output_file, str(out_path))

    def test_raises_error_if_no_cram_index_path_given(self, tmp_path: Path):
        _, cram_pth_no_index, batch = setup_test(tmp_path, index=False)

        with pytest.raises(ValueError, match='CRAM for somalier is required to have CRAI index'):
            _ = extract(
                b=batch,
                cram_path=cram_pth_no_index,
                out_somalier_path=(tmp_path / 'output_file'),
            )
