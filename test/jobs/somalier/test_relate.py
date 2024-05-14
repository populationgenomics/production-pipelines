import re

import pytest
from pytest_mock import MockFixture

from cpg_utils import Path
from cpg_workflows.jobs.somalier import MAX_FREEMIX, pedigree

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
    dataset.add_sequencing_group(id='CPGAAAAAA', external_id='SAMPLE1')
    batch = create_local_batch(tmp_path)

    sg_id = dataset.get_sequencing_group_ids()[0]
    somalier_path_by_sgid = {sg_id: (tmp_path / 'test.somalier')}

    return config, batch, somalier_path_by_sgid, dataset


class TestSomalierRelate:
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

        relate_jobs = batch.select_jobs('Somalier relate')
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
        relate_j = pedigree_jobs[0]

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
        relate_j = pedigree_jobs[0]

        assert relate_j is not None
        assert relate_j.attributes == expected_attrs

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
        relate_j = pedigree_jobs[0]

        assert relate_j is not None
        assert relate_j._image == config.images['somalier']

    def test_if_verifybamid_exists_for_sg_check_freemix(self, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

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

        cmd = get_command_str(relate_j)
        sg = dataset.get_sequencing_groups()[0]
        verifybamid_file = somalier_path_by_sgid[sg.id].name  # same as somalier file
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

    def test_if_verifybamid_file_not_available_inputs_files_still_generated(self, tmp_path: Path):
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
        relate_j = pedigree_jobs[0]

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
        relate_j = pedigree_jobs[0]

        cmd = get_command_str(relate_j)
        expected_ped_filename = 'test_ped.ped'
        sample_id_list_file = 'sample_ids.list'
        assert re.search(fr'cat \${{BATCH_TMPDIR}}/inputs/\w+/{expected_ped_filename}', cmd)
        assert re.search(fr'grep -f \$BATCH_TMPDIR/{sample_id_list_file}', cmd)

    def test_moves_somalier_output_to_expected_batch_resource_files(self, tmp_path: Path):
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
        relate_j = pedigree_jobs[0]

        cmd = get_command_str(relate_j)
        assert re.search(
            r'mv related.pairs.tsv \${BATCH_TMPDIR}/.+\/output_pairs',
            cmd,
        )
        assert re.search(
            r'mv related.samples.tsv \${BATCH_TMPDIR}/.+\/output_samples',
            cmd,
        )
        assert re.search(
            r'mv related.html \${BATCH_TMPDIR}/.+\/output_html',
            cmd,
        )

    def test_related_html_sequencing_group_id_replacement(self, tmp_path: Path):
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
        relate_j = pedigree_jobs[0]

        cmd = get_command_str(relate_j)
        sequencing_group_id = dataset.get_sequencing_groups()[0].id
        rich_id = dataset.rich_id_map()[sequencing_group_id]

        assert re.search(
            fr"sed -iBAK 's/{sequencing_group_id}/{rich_id}/g' related.html",
            cmd,
        )

    def test_writes_outputs_to_final_destinations(self, mocker: MockFixture, tmp_path: Path):
        _, batch, somalier_path_by_sgid, dataset = setup_pedigree_test(tmp_path)

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

        assert relate_j is not None
        spy.assert_has_calls(
            calls=[
                mocker.call(relate_j.output_samples, str(out_samples_path)),
                mocker.call(relate_j.output_pairs, str(out_pairs_path)),
                mocker.call(relate_j.output_html, str(out_html_path)),
            ],
        )
