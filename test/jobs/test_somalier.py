import re

import pytest
from cpg_utils import Path
from pytest_mock import MockFixture

from cpg_workflows.jobs.somalier import MAX_FREEMIX, extract, pedigree

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.dataset import create_dataset
from .helpers import get_command_str


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


def setup_pedigree_test(tmp_path: Path, config: PipelineConfig | None = None):
    config = config or default_config()
    set_config(config, tmp_path / 'config.toml')

    dataset_id = config.workflow.dataset

    dataset = create_dataset(name=dataset_id)
    dataset.add_sequencing_group(id='CPG000001', external_id='SAMPLE1')
    batch = create_local_batch(tmp_path)

    sg_id = dataset.get_sequencing_group_ids()[0]
    somalier_path_by_sgid = {sg_id: (tmp_path / 'test.somalier')}

    return config, batch, somalier_path_by_sgid, dataset


class TestSomalierExtract:
    def test_creates_one_job(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        _ = extract(b=batch, cram_path=cram_pth, out_somalier_path=tmp_path)

        # ---- Assertions
        assert (
            len(batch.select_jobs('Somalier extract')) == 1
        ), "Unexpected number of 'Somalier' jobs in batch list, should be just 1 job"

    def test_will_return_none_if_output_file_exists_and_overwrite_is_false(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=False,
        )

        # --- Assertions
        assert j is None, 'Job was created when output file should have been reused'

    def test_will_create_one_job_if_output_file_exists_and_overwrite_is_true(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        output = tmp_path / 'output_somalier_file'
        output.touch()

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=tmp_path,
            overwrite=True,
        )

        # ---- Assertions
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
        self, tmp_path: Path, job_attrs, expected_attrs
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
            job_attrs=job_attrs,
        )

        # ---- Assertions
        assert j is not None
        assert j.attributes == expected_attrs

    def test_uses_image_specified_in_config(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        assert j is not None
        assert j._image == config.images['somalier']

    def test_sets_sites_location_with_name_of_sites_file_in_config(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        sites = config.other['references']['somalier_sites']
        assert j is not None
        assert re.search(fr'SITES=\$BATCH_TMPDIR/sites/{sites}', cmd)

    def test_runs_with_somalier_sites_from_config_file(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        sites = config.other['references']['somalier_sites']
        assert j is not None
        assert re.search(fr'--sites \${{BATCH_TMPDIR}}/inputs/\w+/{sites}', cmd)

    def test_uses_fail_safe_copy_on_cram_path_and_index_in_bash_command(
        self, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        assert re.search(fr'retry_gs_cp .*{cram_pth.path}', cmd)
        assert re.search(fr'retry_gs_cp .*{cram_pth.index_path}', cmd)

    def test_uses_reference_in_workflow_config_section_if_set(self, tmp_path: Path):
        # ---- Test setup
        config = default_config()
        config.workflow.ref_fasta = 'test_workflow_ref.fa'
        _, cram_pth, batch = setup_test(tmp_path, config)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        ref_file = config.workflow.ref_fasta
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_uses_broad_reference_as_default_if_reference_not_set_in_workflow_config_section(
        self, tmp_path: Path
    ):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=(tmp_path / 'output_file'),
        )

        # ---- Assertions
        cmd = get_command_str(j)
        ref_file = config.other['references']['broad']['ref_fasta']
        assert re.search(fr'-f \${{BATCH_TMPDIR}}/inputs/\w+/{ref_file}', cmd)

    def test_batch_writes_output_file_to_output_path(
        self, mocker: MockFixture, tmp_path: Path
    ):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path)

        # ---- The jobs we want to test
        spy = mocker.spy(batch, 'write_output')
        out_path = tmp_path / 'output_file'

        j = extract(
            b=batch,
            cram_path=cram_pth,
            out_somalier_path=out_path,
        )

        # ---- Assertions
        assert j is not None
        spy.assert_called_with(j.output_file, str(out_path))

    def test_raises_error_if_no_cram_index_path_given(self, tmp_path: Path):
        # ---- Test setup
        _, cram_pth, batch = setup_test(tmp_path, index=False)

        # ---- The jobs we want to test and Assertions
        with pytest.raises(
            ValueError, match='CRAM for somalier is required to have CRAI index'
        ):
            _ = extract(
                b=batch,
                cram_path=cram_pth,
                out_somalier_path=(tmp_path / 'output_file'),
            )


class TestSomalierPedigree:
    class TestSomalierRelate:
        def test_creates_one_relate_job(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)
            # ---- The job that we want to test
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

            # ---- Assertions
            relate_jobs = batch.select_jobs('Somalier relate')
            assert len(relate_jobs) == 1

        def test_if_label_provided_adds_label_to_default_job_title_(
            self, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            label = 'test-label'

            # ---- The job that we want to test
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
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            job_name = f'Somalier relate [{label}]'
            assert relate_j.name == job_name

        @pytest.mark.parametrize(
            'job_attrs, expected_attrs',
            [
                (None, {'tool': 'somalier'}),
                (
                    {'test_tool': 'test_relate'},
                    {'test_tool': 'test_relate', 'tool': 'somalier'},
                ),
            ],
        )
        def test_if_job_attrs_supplied_job_attrs_set_or_sets_default_attrs_if_not_supplied(
            self, tmp_path: Path, job_attrs, expected_attrs
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            assert relate_j is not None
            assert relate_j.attributes == expected_attrs

        def test_uses_image_specified_in_config(self, tmp_path: Path):
            # ---- Test setup
            config, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(
                tmp_path
            )

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            assert relate_j is not None
            assert relate_j._image == config.images['somalier']

        def test_if_verifybamid_exists_for_sg_check_freemix(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
                verifybamid_by_sgid=somalier_path_by_sgid,
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            cmd = get_command_str(relate_j)
            sg = dataset.get_sequencing_groups()[0]
            verifybamid_file = somalier_path_by_sgid[
                sg.id
            ].name  # same as somalier file
            assert re.search(
                fr'FREEMIX=\$\(cat \${{BATCH_TMPDIR}}/inputs/\w+/{verifybamid_file}',
                cmd,
            )
            assert re.search(fr'\(echo "\$FREEMIX > {MAX_FREEMIX}" \| bc\) -eq 0', cmd)

            # Test somalier file and sg id's appended to correct files
            relate_input_file = 'input_files.list'
            somalier_file = somalier_path_by_sgid[sg.id].name
            sample_id_list_file = 'sample_ids.list'
            sequencing_group_id = sg.id
            assert re.search(
                fr'echo "\${{BATCH_TMPDIR}}/inputs/\w+/{somalier_file}" >> \$BATCH_TMPDIR/{relate_input_file}',
                cmd,
            )
            assert re.search(
                fr'echo "{sequencing_group_id}" >> \$BATCH_TMPDIR/{sample_id_list_file}',
                cmd,
            )

        def test_if_verifybamid_file_not_available_somalier_file_written_to_relate_input_file_(
            self, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            cmd = get_command_str(relate_j)
            sg = dataset.get_sequencing_groups()[0]
            relate_input_file = 'input_files.list'
            somalier_file = somalier_path_by_sgid[sg.id].name
            sample_id_list_file = 'sample_ids.list'
            sequencing_group_id = sg.id
            assert 'FREEMIX' not in cmd
            assert re.search(
                fr'echo "\${{BATCH_TMPDIR}}/inputs/\w+/{somalier_file}" >> \$BATCH_TMPDIR/{relate_input_file}',
                cmd,
            )
            assert re.search(
                fr'echo "{sequencing_group_id}" >> \$BATCH_TMPDIR/{sample_id_list_file}',
                cmd,
            )

        def test_filter_expected_pedigrees_to_ped_file(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            cmd = get_command_str(relate_j)
            expected_ped_filename = 'test_ped.ped'
            sample_id_list_file = 'sample_ids.list'
            assert re.search(
                fr'cat \${{BATCH_TMPDIR}}/inputs/\w+/{expected_ped_filename}', cmd
            )
            assert re.search(fr'grep -f \$BATCH_TMPDIR/{sample_id_list_file}', cmd)

        def test_moves_somalier_output_to_expected_batch_resource_files(
            self, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            cmd = get_command_str(relate_j)
            assert re.search(
                fr'mv related.pairs.tsv \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_pairs',
                cmd,
            )
            assert re.search(
                fr'mv related.samples.tsv \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_samples',
                cmd,
            )
            assert re.search(
                fr'mv related.html \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_html',
                cmd,
            )

        def test_related_html_sequencing_group_id_replacement(
            self, mocker: MockFixture, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            cmd = get_command_str(relate_j)
            sequencing_group_id = dataset.get_sequencing_groups()[0].id
            rich_id = dataset.rich_id_map()[sequencing_group_id]

            assert re.search(
                fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' related.html",
                cmd,
            )

        def test_writes_outputs_to_final_destinations(
            self, mocker: MockFixture, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            spy = mocker.spy(batch, 'write_output')
            out_samples_path = tmp_path / 'out_samples'
            out_pairs_path = tmp_path / 'out_pairs'
            out_html_path = tmp_path / 'out_html'

            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=out_samples_path,
                out_pairs_path=out_pairs_path,
                out_html_path=out_html_path,
            )
            relate_j = pedigree_jobs[0]

            # ---- Assertions
            assert relate_j is not None
            spy.assert_has_calls(
                calls=[
                    mocker.call(relate_j.output_samples, str(out_samples_path)),
                    mocker.call(relate_j.output_pairs, str(out_pairs_path)),
                    mocker.call(relate_j.output_html, str(out_html_path)),
                ]
            )

    class TestSomalierCheckPedigree:
        # Test _check_pedigree
        def test_creates_one_relate_job(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)
            # ---- The job that we want to test
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

            # ---- Assertions
            relate_jobs = batch.select_jobs('Pedigree check')
            assert len(relate_jobs) == 1

        def test_if_label_provided_adds_label_to_default_job_title_(
            self, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            label = 'test-label'

            # ---- The job that we want to test
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

            # ---- Assertions
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
            self, tmp_path: Path, job_attrs, expected_attrs
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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

            # ---- Assertions
            assert pedigree_check_j is not None
            assert pedigree_check_j.attributes == expected_attrs

        def test_uses_image_specified_in_config(self, tmp_path: Path):
            # ---- Test setup
            config, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(
                tmp_path
            )

            # ---- The job that we want to test
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

            # ---- Assertions
            assert check_pedigree_j is not None
            assert check_pedigree_j._image == config.images['cpg_workflows']

        def test_script_path_correctness(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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
            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            script_name = 'check_pedigree.py'
            assert re.search(fr'python3 {script_name}', cmd)

        def test_sed_commands_for_samples_pairs_and_expected_ped_files(
            self, tmp_path: Path
        ):
            # ---- Test setup
            config, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(
                tmp_path
            )

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]
            check_pedigree_j = pedigree_jobs[1]

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            sequencing_group_id = dataset.get_sequencing_groups()[0].id
            rich_id = dataset.rich_id_map()[sequencing_group_id]
            expected_ped = 'test_ped.ped'

            # check_pedigree uses output_samples and output_pairs files from _relate job
            assert re.search(
                fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_samples'",
                cmd,
            )
            assert re.search(
                fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_pairs'",
                cmd,
            )
            assert re.search(
                fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' \${{BATCH_TMPDIR}}/inputs/\w+/{expected_ped}",
                cmd,
            )

        def test_sed_commands_not_present_if_no_rich_id_map_provided(
            self, tmp_path: Path
        ):
            # ---- Test setup
            config, batch, somalier_path_by_sgid, _ = setup_pedigree_test(tmp_path)

            dataset_id = config.workflow.dataset
            dataset = create_dataset(name=dataset_id)
            dataset.add_sequencing_group(id='CPG000001')
            # ---- The job that we want to test
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

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)

            assert 'sed -iBAK' not in cmd
            assert bool(dataset.rich_id_map()) is False  # empty dictionary

        def test_flags_called_correctly(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
            )
            relate_j = pedigree_jobs[0]
            check_pedigree_j = pedigree_jobs[1]

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            expected_ped = 'test_ped.ped'
            somalier_html_url = 'test_html_url'
            assert re.search(
                fr'--somalier-samples \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_samples',
                cmd,
            )
            assert re.search(
                fr'--somalier-pairs \${{BATCH_TMPDIR}}/{relate_j._dirname}/output_pairs',
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
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            assert re.search('--no-send-to-slack', cmd)

        def test_if_out_html_url_flag_is_set(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
            pedigree_jobs = pedigree(
                dataset=dataset,
                b=batch,
                expected_ped_path=(tmp_path / 'test_ped.ped'),
                somalier_path_by_sgid=somalier_path_by_sgid,
                out_samples_path=(tmp_path / 'out_samples'),
                out_pairs_path=(tmp_path / 'out_pairs'),
                out_html_path=(tmp_path / 'out_html'),
                out_html_url='test_html_url',
                send_to_slack=False,
            )
            check_pedigree_j = pedigree_jobs[1]

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            somalier_html_url = 'test_html_url'
            assert re.search(fr'--html-url {somalier_html_url}', cmd)

        def test_job_output_created(self, tmp_path: Path):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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

            # ---- Assertions
            cmd = get_command_str(check_pedigree_j)
            job_output = 'output'
            assert re.search(
                fr'touch \${{BATCH_TMPDIR}}/{check_pedigree_j._dirname}/{job_output}',
                cmd,
            )

        def test_if_output_path_provided_writes_outputs_to_final_destination(
            self, mocker: MockFixture, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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

            # ---- Assertions
            out_checks_path = tmp_path / 'out_path'
            assert check_pedigree_j is not None
            spy.assert_has_calls(
                calls=[mocker.call(check_pedigree_j.output, str(out_checks_path))]
            )

        def test_if_output_path_not_provided_writes_outputs_to_final_destination(
            self, mocker: MockFixture, tmp_path: Path
        ):
            # ---- Test setup
            _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

            # ---- The job that we want to test
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

            # ---- Assertions
            out_checks_path = tmp_path / 'out_path'
            assert check_pedigree_j is not None
            with pytest.raises(AssertionError):
                spy.assert_has_calls(
                    calls=[mocker.call(check_pedigree_j.output, str(out_checks_path))],
                    any_order=True,
                )
