"""
Test the trim stage
"""

from pathlib import Path
from cpg_utils import to_path
from unittest.mock import mock_open
from pytest_mock import MockFixture

from .. import update_dict, set_config

from os import makedirs
import os.path

def get_toml(tmp_path) -> str:
    return f"""
    [workflow]
    name = "rare_rna_trim"
    dataset_gcp_project = "test-analysis-dataset-1234"
    dataset = "test-analysis-dataset"
    access_level = "test"
    sequencing_type = "rna"
    driver_image = "<stub>"
    check_inputs = false
    check_intermediates = false
    check_expected_outputs = false
    path_scheme = "local"
    local_dir = "{tmp_path}"

    [hail]
    billing_project = "test-analysis-dataset"
    delete_scratch_on_exit = true
    dry_run = true
    backend = "local"

    [images]
    cutadapt = "stub"

    [storage.default]
    default = '{tmp_path}'
    web = "{tmp_path}-web"
    analysis = "{tmp_path}-analysis"
    tmp = "{tmp_path}-test-tmp"
    web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

    [storage.test-input-dataset]
    default = "{tmp_path}"
    web = "{tmp_path}-web"
    analysis = "{tmp_path}-analysis"
    tmp = "{tmp_path}-test-tmp"
    web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

    [storage.test-analysis-dataset]
    default = "{tmp_path}"
    web = "{tmp_path}-web"
    analysis = "{tmp_path}-analysis"
    tmp = "{tmp_path}-test-tmp"
    web_url = "https://test-web.populationgenomics.org.au/fewgenomes"
    """


DEFAULT_CONFIG = Path(
    to_path(__file__).parent.parent.parent / 'cpg_workflows' / 'defaults.toml'
)


def _mock_cohort():
    from cpg_workflows.targets import Cohort
    from cpg_workflows.filetypes import FastqPair, FastqPairs  # , BamPath

    cohort = Cohort()
    ds = cohort.create_dataset('test-input-dataset')
    ds.add_sequencing_group(
        'CPG01',
        'SAMPLE1',
        alignment_input_by_seq_type={
            'rna': FastqPairs(
                [
                    FastqPair(
                        'gs://test-input-dataset-upload/sample1_L1_R1.fq.gz',
                        'gs://test-input-dataset-upload/sample1_L1_R2.fq.gz',
                    ),
                    FastqPair(
                        'gs://test-input-dataset-upload/sample1_L2_R1.fq.gz',
                        'gs://test-input-dataset-upload/sample1_L2_R2.fq.gz',
                    ),
                ]
            )
        },
    )
    return cohort


def selective_mock_open(*args, **kwargs):
    if str(args[0]).endswith('.toml'):
        # Don't mock calls to load a config file
        return open(*args, **kwargs)
    else:
        return mock_open(read_data='<stub>')(*args, **kwargs)


def test_rare_rna(mocker: MockFixture, tmp_path):
    """
    Test the trim stage in a dry run mode.
    """
    from hailtop.batch.job import Job
    from cpg_workflows.jobs.trim import trim

    conf = get_toml(tmp_path)
    set_config(
        conf,
        tmp_path / 'config.toml',
        merge_with=[DEFAULT_CONFIG],
    )

    mocker.patch('cpg_workflows.inputs.create_cohort', _mock_cohort)

    def do_nothing(*args, **kwargs):
        return None

    def mock_create_analysis(*args, **kwargs) -> int:
        return 1
    
    # Capture the trim command
    cmd_str_list = []
    def capture_trim_cmd(*args, **kwargs) -> list[Job]:
        trim_job = trim(*args, **kwargs)
        cmd_str_list.append(
            "===== FASTQ TRIM JOB START =====\n\n" +
            "\n".join(trim_job._command) +
            "\n\n===== FASTQ TRIM JOB END =====\n\n"
        )
        return trim_job

    mocker.patch('pathlib.Path.open', selective_mock_open)
    # functions like get_intervals checks file existence
    mocker.patch('cpg_workflows.workflow.list_all_parent_dirs', lambda *args: {})
    mocker.patch('cpg_workflows.workflow.list_of_all_dir_contents', lambda *args: {})
    mocker.patch('cpg_workflows.workflow.missing_from_pre_collected', lambda *args: None)
    # cloudfuse (used in Vep) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.cloudfuse', do_nothing)
    # always_run (used in MtToEs -> hail_dataproc_job) doesn't work with LocalBackend
    mocker.patch('hailtop.batch.job.Job.always_run', do_nothing)
    # can't access secrets from CI environment
    mocker.patch(
        'cpg_workflows.stages.seqr_loader.es_password', lambda: 'test-password'
    )
    mocker.patch(
        'metamist.apis.AnalysisApi.create_analysis',
        mock_create_analysis,
    )
    mocker.patch('metamist.apis.AnalysisApi.update_analysis', do_nothing)

    # Path the trim function to capture the job command
    mocker.patch('cpg_workflows.jobs.trim.trim', capture_trim_cmd)

    from cpg_workflows.batch import get_batch
    from cpg_workflows.inputs import get_cohort
    from cpg_workflows.workflow import get_workflow

    # Imports specific for testing the trim stage
    from cpg_utils.hail_batch import dataset_path
    from cpg_workflows.stages.trim import Trim
    from cpg_workflows.filetypes import FastqPairs, FastqPair

    get_workflow().run(stages=[Trim])

    b = get_batch()
    trim_job = b.job_by_tool['cutadapt']
    sample_list = get_cohort().get_sequencing_groups()

    # The number of FASTQ trim jobs should equal the number of FASTQ pairs
    n_trim_jobs_list = [
        len(s.alignment_input_by_seq_type.get('rna'))
        for s in sample_list
        if isinstance(s.alignment_input_by_seq_type.get('rna'), FastqPairs)
    ]
    n_trim_jobs_list.extend([
        1
        for s in sample_list
        if isinstance(s.alignment_input_by_seq_type.get('rna'), FastqPair)
    ])
    n_trim_jobs = sum(n_trim_jobs_list)
    assert len(trim_job) == n_trim_jobs

    output_path = dataset_path('cmd.txt')
    output_path_parent = os.path.dirname(output_path)
    makedirs(output_path_parent, exist_ok=True)
    with open(output_path, 'w') as f:
        for cmd_str in cmd_str_list:
            f.write(cmd_str)
        f.write('\n')
