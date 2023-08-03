from cpg_utils import Path
from cpg_utils.hail_batch import image_path
from pytest_mock import MockFixture

from cpg_workflows.jobs.verifybamid import verifybamid

from .. import set_config
from ..factories.alignment_input import create_cram_input
from ..factories.batch import create_local_batch
from ..factories.config import PipelineConfig, WorkflowConfig
from ..factories.sequencing_group import create_sequencing_group
from .helpers import get_command_str


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='verifybamid-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'verifybamid': 'test_image',
        },
        other={
            'references': {
                'broad': {
                    'ref_fasta': 'hg38_reference.fa',
                    'dragmap_prefix': 'gs://a-cpg-bucket/dragen_reference/',
                    'genome_contam_ud': 'test_genome_ud.ud',
                    'genome_contam_bed': 'test_genome_bed.bed',
                    'genome_contam_mu': 'test__genome_mu.mu',
                    'exome_contam_ud': 'test_exome_ud.ud',
                    'exome_contam_bed': 'test_exome_bed.bed',
                    'exome_contam_mu': 'test_exome_mu.mu',
                }
            },
            'cramqc': {'num_pcs': '4'},
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


class TestVerifyBAMID:
    def test_VBI_creates_a_job(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=tmp_path,
            job_attrs=None,
            overwrite=True,
        )

        # ---- Assertions
        assert (
            j is not None
        ), 'The verifybamid function did not create a job. Check if overwrite=False in verifybamid call'

    def test_VBI_can_reuse_existing_output_path(self, tmp_path: Path):
        # ---- Test setup
        config, cram_pth, batch = setup_test(tmp_path)

        # ---- The job we want to test
        j = verifybamid(
            b=batch,
            cram_path=cram_pth,
            out_verify_bamid_path=tmp_path,
            job_attrs=None,
            overwrite=False,
        )
        print()

    def test_VBI_can_overwrite_existing_output_path(self, tmp_path: Path):
        pass

    def test_sets_job_attrs_or_sets_default_attrs_if_not_supplied(self, tmp_path: Path):
        pass

    def test_uses_VBI_image_specified_in_config(self, tmp_path: Path):
        pass
