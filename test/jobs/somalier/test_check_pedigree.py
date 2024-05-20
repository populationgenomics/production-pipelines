import re

import pytest
from pytest_mock import MockFixture

from cpg_utils import Path
from cpg_workflows.jobs.somalier import pedigree

from ... import set_config
from ...factories.batch import create_local_batch
from ...factories.config import PipelineConfig, WorkflowConfig
from ...factories.dataset import create_dataset
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


def setup_pedigree_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    dataset_id = config.workflow.dataset

    dataset = create_dataset(name=dataset_id)
    dataset.add_sequencing_group(id='CPGABCD', external_id='SAMPLE1')
    batch = create_local_batch(tmp_path)

    sg_id = dataset.get_sequencing_group_ids()[0]
    somalier_path_by_sgid = {sg_id: (tmp_path / 'test.somalier')}

    return config, batch, somalier_path_by_sgid, dataset


class TestSomalierCheckPedigree:
    def test_creates_one_relate_job(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        _ = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            label=None,
        )

        relate_jobs = batch.select_jobs('Pedigree check')
        assert len(relate_jobs) == 1

    def test_if_label_provided_adds_label_to_default_job_title_(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        label = 'test-label'

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            label=label,
        )
        pedigree_check_j = pedigree_jobs[1]

        job_name = f'Pedigree check [{label}]'
        assert pedigree_check_j.name == job_name

    @pytest.mark.parametrize(
        'job_attrs, expected_attrs',
        [
            (None, {'tool': 'python'}),
            (
                {'test_tool': 'test_check_pedigree'},
                {'test_tool': 'test_check_pedigree', 'tool': 'python'},
            ),
        ],
    )
    def test_if_job_attrs_supplied_job_attrs_set_or_sets_default_attrs_if_not_supplied(
        self,
        tmp_path: Path,
        job_attrs,
        expected_attrs,
    ):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            job_attrs=job_attrs,
        )
        pedigree_check_j = pedigree_jobs[1]

        assert pedigree_check_j is not None
        assert pedigree_check_j.attributes == expected_attrs

    def test_uses_image_specified_in_config(self, tmp_path: Path):
        config, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]

        assert check_pedigree_j is not None
        assert check_pedigree_j._image == config.images['cpg_workflows']

    def test_script_path_correctness(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]
        cmd = get_command_str(check_pedigree_j)
        script_name = 'check_pedigree.py'
        assert re.search(fr'python3 {script_name}', cmd)

    def test_sed_commands_for_samples_pairs_and_expected_ped_files(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)
        sequencing_group_id = dataset.get_sequencing_groups()[0].id
        rich_id = dataset.rich_id_map()[sequencing_group_id]
        expected_ped = 'test_ped.ped'

        assert re.search(
            fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/.+\/output_samples'",
            cmd,
        )
        assert re.search(
            fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/.+\/output_pairs'",
            cmd,
        )
        assert re.search(
            fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/inputs/\w+/{expected_ped}",
            cmd,
        )

    def test_sed_commands_not_present_if_no_rich_id_map_provided(self, tmp_path: Path):
        config, batch, somalier_path_by_sgid, _ = setup_pedigree_test(tmp_path)

        dataset_id = config.workflow.dataset
        dataset = create_dataset(name=dataset_id)
        dataset.add_sequencing_group(id='CPGABCD')

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)

        assert 'sed -iBAK' not in cmd
        assert bool(dataset.rich_id_map()) is False  # empty dictionary

    def test_check_pedigree_uses_output_files_from_relate(self, mocker: MockFixture, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        relate_j, check_pedigree_j = pedigree_jobs
        relate_cmd = get_command_str(relate_j)
        check_pedigree_j_cmd = get_command_str(check_pedigree_j)
        out_sample_path_match = re.search(r'\${BATCH_TMPDIR}\/.+\/output_samples', relate_cmd)
        out_pairs_path_match = re.search(r'\${BATCH_TMPDIR}\/.+\/output_pairs', relate_cmd)
        if out_sample_path_match is not None and out_pairs_path_match is not None:
            matched_sample = out_sample_path_match.group(0)
            matched_pairs = out_pairs_path_match.group(0)

        assert f'--somalier-samples {matched_sample}' in check_pedigree_j_cmd
        assert f'--somalier-pairs {matched_pairs}' in check_pedigree_j_cmd

    def test_flags_called_correctly(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)
        expected_ped = 'test_ped.ped'
        somalier_html_url = 'test_html_url'
        assert re.search(
            r'--somalier-samples \${BATCH_TMPDIR}\/.+\/output_samples',
            cmd,
        )
        assert re.search(
            r'--somalier-pairs \${BATCH_TMPDIR}\/.+\/output_pairs',
            cmd,
        )
        assert re.search(
            fr'--expected-ped \${{BATCH_TMPDIR}}/inputs/\w+/{expected_ped}',
            cmd,
        )
        assert not re.search(fr'--html-url {somalier_html_url}', cmd)
        assert re.search(fr'--dataset {dataset.name}', cmd)
        assert re.search(fr'--title "{check_pedigree_j.name}"', cmd)
        assert re.search('--send-to-slack', cmd)

    def test_if_send_to_slack_false_flag(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            send_to_slack=False,
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)
        assert re.search('--no-send-to-slack', cmd)

    @pytest.mark.parametrize('out_html_url', ['test_html_url', None])
    def test_if_out_html_url_flag_is_set(self, tmp_path: Path, out_html_url: str):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            out_html_url=out_html_url,
            send_to_slack=False,
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)
        pattern = r'python3 check_pedigree\.py.*?--html-url ([^\n]+)'
        regex = re.compile(pattern, re.DOTALL)
        match = regex.search(cmd)
        if out_html_url is None:
            assert bool(match) is False
        else:
            assert re.search(fr'--html-url {out_html_url}', cmd)

    def test_job_output_created(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
        )
        check_pedigree_j = pedigree_jobs[1]

        cmd = get_command_str(check_pedigree_j)
        job_output = 'output'
        assert re.search(
            fr'touch \${{BATCH_TMPDIR}}/{check_pedigree_j._dirname}/{job_output}',
            cmd,
        )

    def test_if_output_path_provided_writes_outputs_to_final_destination(self, mocker: MockFixture, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        spy = mocker.spy(batch, 'write_output')

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            out_checks_path=(tmp_path / 'out_path'),
        )
        check_pedigree_j = pedigree_jobs[1]

        out_checks_path = tmp_path / 'out_path'
        assert check_pedigree_j is not None
        spy.assert_has_calls(calls=[mocker.call(check_pedigree_j.output, str(out_checks_path))])

    def test_if_output_path_not_provided_does_not_write_outputs_to_final_destination(
        self,
        mocker: MockFixture,
        tmp_path: Path,
    ):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

        spy = mocker.spy(batch, 'write_output')

        pedigree_jobs = pedigree(
            dataset=dataset,
            b=batch,
            expected_ped_path=(tmp_path / 'test_ped.ped'),
            somalier_path_by_sgid=somalier_path_by_sgid,
            out_samples_path=(tmp_path / 'out_samples'),
            out_pairs_path=(tmp_path / 'out_pairs'),
            out_html_path=(tmp_path / 'out_html'),
            # No output path provided
        )
        check_pedigree_j = pedigree_jobs[1]

        out_checks_path = tmp_path / 'out_path'
        assert check_pedigree_j is not None
        with pytest.raises(AssertionError):
            spy.assert_has_calls(
                calls=[mocker.call(check_pedigree_j.output, str(out_checks_path))],
                any_order=True,
            )
