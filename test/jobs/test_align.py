"""
Testing Alignment Job
"""

from cpg_utils import to_path

from .. import set_config

from pytest_mock import MockFixture

from cpg_workflows.batch import get_batch
from cpg_workflows.inputs import _populate_alignment_inputs
from cpg_workflows.targets import SequencingGroup, Dataset

from cpg_workflows.jobs.align import align

example_entry = {
    'type': 'genome',
    'assays': [
        {
            'id': 8,
            'meta': {
                'concentration': '25.8',
                'reference_genome': 'Homo sapiens (b37d5)',
                'fluid_x_tube_id': '220405_FS28686769',
                'volume': '100',
                'platform': '30x Illumina PCR-Free',
                'reads_type': 'fastq',
                'batch': 'M001',
                'reads': [
                    {
                        'location': 'gs://cpg-fewgenomes-main-upload/EXAMPLE_M001_R1.fastq.gz',
                        'basename': 'EXAMLPE_M001_R1.fastq.gz',
                        'class': 'File',
                        'checksum': 'md5:dbb1d8db1e08695a9d6b1212974a6bd7  EXAMPLE_M001_R1.fastq.gz',
                        'size': 8171768680,
                        'datetime_added': None,
                    },
                    {
                        'location': 'gs://cpg-fewgenomes-main-upload/EXAMPLE_M001_R2.fastq.gz',
                        'basename': 'EXAMPLE_M001_R2.fastq.gz',
                        'class': 'File',
                        'checksum': 'md5:937fd847029f4f5ea254df397efd1c31  EXAMPLE_M001_R2.fastq.gz',
                        'size': 8773022465,
                        'datetime_added': None,
                    },
                ],
                'sequencing_type': 'genome',
                'sequencing_technology': 'short-read',
                'sequencing_platform': 'illumina',
            },
            'type': 'sequencing',
        }
    ],
}


def test_align_job(mocker: MockFixture, tmp_path):
    config = f"""
    [workflow]
    dataset_gcp_project = 'fewgenomes'
    access_level = 'test'
    dataset = 'fewgenomes'
    driver_image = '<stub>'
    sequencing_type = 'genome'
    status_reporter = 'metamist'

    check_inputs = false
    check_intermediates = false
    check_expected_outputs = false
    path_scheme = 'local'

    [storage.default]
    default = '{tmp_path}'

    [storage.fewgenomes]
    default = '{tmp_path}'

    [hail]
    billing_project = 'fewgenomes'
    delete_scratch_on_exit = false
    backend = 'local'
    dry_run = true

    [resource_overrides]
    # Override default resource requirements for unusually large seq data without
    # demanding higher resources for all operations as standard. Examples below


    [images]
    dragmap = 'australia-southeast1-docker.pkg.dev/cpg-common/images/dragmap:1.3.0'
    picard = 'australia-southeast1-docker.pkg.dev/cpg-common/images/picard:2.27.4'

    [references.broad]
    dragmap_prefix = 'gs://cpg-common-main/references/hg38/v0/dragen_reference'
    ref_fasta = 'gs://cpg-common-main/references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta'
    """

    set_config(config, tmp_path / 'config.toml')

    fewgenomes_dataset = Dataset(name='fewgenomes')

    sequencing_group = SequencingGroup(
        id='CPG01', external_id='SAMPLE1', dataset=fewgenomes_dataset
    )

    # Is this bad to do? It feels bad. Should we just generate the cohort?
    _populate_alignment_inputs(sequencing_group=sequencing_group, entry=example_entry)

    align_jobs = align(b=get_batch(), sequencing_group=sequencing_group)

    assert len(align_jobs) == 2

    for job in align_jobs:
        print(job._command)
        # TODO: Best way to test this??
