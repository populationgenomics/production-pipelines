import re

import pytest
from cpg_utils import Path
from cpg_utils.config import ConfigError
from pytest_mock import MockFixture

from cpg_workflows.jobs import multiqc as mqc_module
from cpg_workflows.jobs.multiqc import (
    _write_sg_id_map,
    check_multiqc,
    check_report_job,
    multiqc,
)

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.dataset import create_dataset
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='multiqc-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
            dry_run=False,  # may need to change this
        ),
        images={
            'multiqc': 'test_image',
            'cpg_workflows': 'test_check_report_job_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                },
            },
            'qc_thresholds': {
                'genome': {'min': {'genome_metric1': 20, 'genome_metric2': 20}},
                'exome': {'max': {'exome_metric1': 20, 'exome_metric2': 20}},
            },
        },
    )


def setup_multiqc_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    dataset_id = config.workflow.dataset

    dataset = create_dataset(name=dataset_id)
    dataset.add_sequencing_group(id='CPG000000', external_id='SAMPLE0')
    dataset.add_sequencing_group(id='CPG000001', external_id='SAMPLE1')
    dataset.add_sequencing_group(id='CPG000002', external_id='SAMPLE2')
    dataset.add_sequencing_group(id='CPG000003', external_id='SAMPLE3')
    batch = create_local_batch(tmp_path)

    # Don't think I actually need the following, the endings aren't getting cut as the script isn't running
    paths = [(tmp_path / f'mqc_input_file_CPG00000{i}') for i in range(4)]
    return config, batch, dataset, paths


# def setup_check_report_job_test(tmp_path: Path, config: PipelineConfig | None = None):
#     config = config or default_config()

#     # TODO: put qc_thresholds and cpg_workflows inside default config.
#     if config.images['cpg_workflows'] is not None:
#         config.images['cpg_workflows'] = 'test_check_report_job_image'
#     set_config(config, tmp_path / 'config.toml')

#     dataset_id = config.workflow.dataset
#     dataset = create_dataset(name=dataset_id)
#     dataset.add_sequencing_group(id='CPG000000', external_id='SAMPLE0')
#     dataset.add_sequencing_group(id='CPG000001', external_id='SAMPLE1')
#     dataset.add_sequencing_group(id='CPG000002', external_id='SAMPLE2')
#     dataset.add_sequencing_group(id='CPG000003', external_id='SAMPLE3')

#     batch = create_local_batch(tmp_path)
#     # Don't think I actually need the following, the endings aren't getting cut as the script isn't running
#     paths = [(tmp_path / f'mqc_input_file_CPG00000{i}') for i in range(4)]
#     return config, batch, dataset, paths


class TestMultiQC:
    def test_creates_one_multiqc_job(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        _ = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )

        # ---- Assertions
        mqc_j = batch.select_jobs('MultiQC')
        assert len(mqc_j) == 1

    def test_if_qc_thresholds_provided_check_report_job_created(self, tmp_path: Path):
        # ---- Test setup
        # config = default_config()
        # qc_thresholds = {
        #     'qc_thresholds': {
        #         'genome': {'min': {'genome_metric1': 20, 'genome_metric2': 20}},
        #         'exome': {'max': {'exome_metric1': 20, 'exome_metric2': 20}},
        #     }
        # }
        # config.other['qc_thresholds'] = qc_thresholds
        # config.images['cpg_workflows'] = 'test_check_report_job_image'
        config, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        _ = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )

        # ---- Assertions
        check_j = batch.select_jobs('MultiQC check')
        assert len(check_j) == 1

    def test_if_label_provided_name_of_job_changes(self, tmp_path: Path):
        # ---- Test setup
        config, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        label = 'test-label'
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            label=label,
        )
        mqc_j = jobs[0]
        # ---- Assertions
        mqc_j_name = f'MultiQC [{label}]'
        assert mqc_j.name == mqc_j_name

    def test_mqc_file_list_created(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.dry_run = False
        _, batch, dataset, paths = setup_multiqc_test(tmp_path, config)

        # ---- The job that we want to test
        _ = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )

        # ---- Assertions
        assert (tmp_path / 'multiqc-file-list.txt').exists()

    def test_paths_written_to_mqc_file_list(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.dry_run = False
        _, batch, dataset, paths = setup_multiqc_test(tmp_path, config)

        # ---- The job that we want to test
        _ = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )

        # ---- Assertions
        with open((tmp_path / 'multiqc-file-list.txt'), 'r') as f:
            f_content = f.read()
            for path in paths:
                assert str(path) in f_content

    def test_if_ending_to_trim_provided_file_endings_are_trimmed(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            ending_to_trim=['CPG000000', 'CPG000001', 'CPG000002', 'CPG000003'],
        )
        mqc_j = jobs[0]
        # ---- Assertions
        cmd = get_command_str(mqc_j)
        string_endings_cut = '[CPG000000, CPG000001, CPG000002, CPG000003]'
        assert re.search(
            re.escape(fr'--cl-config "extra_fn_clean_exts: {string_endings_cut}"'), cmd
        )

    def test_if_modules_to_trim_endings_provided_modules_from_string_removed(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            modules_to_trim_endings=['module_0', 'module_1', 'module_2', 'module_3'],
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        modules_to_trim = '[module_0, module_1, module_2, module_3]'
        assert re.search(
            re.escape(fr'--cl-config "use_filename_as_sample_name: {modules_to_trim}"'),
            cmd,
        )

    def test_if_sequencing_group_id_map_provided_sample_map_file_created(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        config = default_config()
        config.workflow.dry_run = False
        _, batch, dataset, paths = setup_multiqc_test(tmp_path, config)
        # TODO: Ask: can i parametrize this if i need access to dataset
        sequencing_group_id_map = {
            sg.target_id: sg.participant_id
            for sg in dataset._sequencing_group_by_id.values()
        }
        sample_map_path = tmp_path / 'rename-sample-map.tsv'
        spy = mocker.spy(mqc_module, '_write_sg_id_map')

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            sequencing_group_id_map=sequencing_group_id_map,
        )
        _ = jobs[0]

        # ---- Assertions
        assert sample_map_path.exists()
        spy.assert_called_once_with(sequencing_group_id_map, sample_map_path)

    def test_if_no_sequencing_group_id_map_provided_sample_map_file_not_created(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config = default_config()
        config.workflow.dry_run = False
        _, batch, dataset, paths = setup_multiqc_test(tmp_path, config)

        sample_map_path = tmp_path / 'rename-sample-map.tsv'

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            sequencing_group_id_map=None,
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        assert not sample_map_path.exists()
        assert '--replace-names' not in cmd

    @pytest.mark.parametrize(
        'extra_config', [{'test_extra1': 'config1', 'test_extra2': 'config2'}, None]
    )
    def test_if_extra_configs_provided_set_as_flags(self, tmp_path: Path, extra_config):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            extra_config=extra_config,
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        serialised_extra_configs = (
            ', '.join(f'{k}: {v}' for k, v in extra_config.items())
            if extra_config
            else ''
        )
        assert (serialised_extra_configs == '') or (
            f'--cl-config "{serialised_extra_configs}"' in cmd
        )

    def test_title_html_and_dataset_names_in_cmd(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        report_filename = 'report'
        dataset_name = dataset.name
        title = mqc_j.name
        assert re.search(fr'--title "{title} for dataset <b>{dataset_name}</b>"', cmd)
        assert re.search(fr'--filename {report_filename}.html', cmd)

    def test_generate_and_copy_multiqc_reports(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        report_filename = 'report'
        assert f'ls output/{report_filename}_data' in cmd
        assert re.search(
            fr'cp output/{report_filename}.html \${{BATCH_TMPDIR}}/\w+-\w+/html', cmd
        )
        assert re.search(
            fr'cp output/{report_filename}_data/multiqc_data.json \${{BATCH_TMPDIR}}/\w+-\w+/json',
            cmd,
        )

    def test_if_out_html_provided_url_written_to_cmd(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)
        out_html_url = 'test_out_html'

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            out_html_url=out_html_url,
        )
        mqc_j = jobs[0]

        # ---- Assertions
        cmd = get_command_str(mqc_j)
        assert f'echo "HTML URL: {out_html_url}' in cmd

    def test_batch_writes_output_to_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)
        out_json_path = tmp_path / 'out_json_path'
        out_html_path = tmp_path / 'out_html_path'
        spy = mocker.spy(batch, 'write_output')

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=out_json_path,
            out_html_path=out_html_path,
        )
        mqc_j = jobs[0]

        # ---- Assertions
        spy.assert_has_calls(
            calls=[
                mocker.call(mqc_j.html, str(out_html_path)),
                mocker.call(mqc_j.json, str(out_json_path)),
            ]
        )

    def test_check_report_job_called_with_correct_params(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)
        spy = mocker.spy(mqc_module, 'check_report_job')

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )
        mqc_j = jobs[0]

        # ---- Assertions
        spy.assert_called_once_with(
            b=batch,
            multiqc_json_file=mqc_j.json,
            multiqc_html_url=None,
            rich_id_map=dataset.rich_id_map(),
            dataset_name=dataset.name,
            label=None,
            out_checks_path=None,
            job_attrs=None,
        )


class TestCheckReport:
    def test_if_label_provided_job_name_changes(self, tmp_path: Path):
        # ---- Test setup
        _, batch, dataset, paths = setup_multiqc_test(tmp_path)
        label = 'test_check_j'

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
            label=label,
        )
        _, check_j = jobs

        # ---- Assertions
        check_j_name = 'MultiQC [test_check_j] check'
        assert check_j.name == check_j_name

    @pytest.mark.parametrize('image', ['test_check_reports_job_image', None])
    def test_if_no_cpg_workflows_image_provided_errors_out(self, tmp_path: Path, image):
        # ---- Test setup
        config = default_config()
        # When writting TOML files, dictionaries with None values are not written, hence why {'cpg_workflows' : None} is not appearing
        config.images['cpg_workflows'] = image
        config, batch, dataset, paths = setup_multiqc_test(tmp_path, config)

        # ---- The job that we want to test
        if image is None:
            with pytest.raises(ConfigError, match='Key "cpg_workflows" not found .*'):
                jobs = multiqc(
                    b=batch,
                    dataset=dataset,
                    tmp_prefix=tmp_path,
                    paths=paths,
                    out_json_path=(tmp_path / 'out_json_path'),
                    out_html_path=(tmp_path / 'out_html_path'),
                )
        if image is not None:
            jobs = multiqc(
                b=batch,
                dataset=dataset,
                tmp_prefix=tmp_path,
                paths=paths,
                out_json_path=(tmp_path / 'out_json_path'),
                out_html_path=(tmp_path / 'out_html_path'),
            )
            _, check_j = jobs
            assert check_j._image == image

    def test_sed_commands_in_cmd(self, tmp_path: Path):
        # ---- Test setup
        config, batch, dataset, paths = setup_multiqc_test(tmp_path)

        # ---- The job that we want to test
        jobs = multiqc(
            b=batch,
            dataset=dataset,
            tmp_prefix=tmp_path,
            paths=paths,
            out_json_path=(tmp_path / 'out_json_path'),
            out_html_path=(tmp_path / 'out_html_path'),
        )
        _, check_j = jobs

        # ---- Assertions
        cmd = get_command_str(check_j)
        rich_id_map = dataset.rich_id_map()
        for sg_id, participant_id in rich_id_map.items():
            assert re.search(fr'sed -iBAK 's/{sg_id}/CPG000000|SAMPLE0/g')
        print()
