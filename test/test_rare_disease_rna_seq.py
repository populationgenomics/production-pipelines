"""
Test the rare disease RNA seq workflow
"""

from pathlib import Path
from cpg_utils import to_path
from unittest.mock import mock_open
from pytest_mock import MockFixture

from . import update_dict, set_config

from os import makedirs
import os.path

import inspect


def get_toml(tmp_path) -> str:
    return f"""
    [workflow]
    name = "rare_disease_rna_seq"
    dataset_gcp_project = "test-analysis-dataset-1234"
    dataset = "test-analysis-dataset"
    access_level = "test"
    sequencing_type = "transcriptome"
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
    fastp = "stub"
    star = "stub"
    samtools = "stub"
    sambamba = "stub"
    subread = "stub"
    outrider = "stub"

    [trim]
    adapter_type = "ILLUMINA_TRUSEQ"

    [references]
    star_ref_dir = "stub"
    gtf = "stub"
    fasta = "stub"

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
    to_path(__file__).parent.parent / 'cpg_workflows' / 'defaults.toml'
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
            'transcriptome': FastqPairs(
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
    ds.add_sequencing_group(
        'CPG02',
        'SAMPLE2',
        alignment_input_by_seq_type={
            'transcriptome': FastqPairs(
                [
                    FastqPair(
                        'gs://test-input-dataset-upload/sample2_L1_R1.fq.gz',
                        'gs://test-input-dataset-upload/sample2_L1_R2.fq.gz',
                    ),
                ]
            ),
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
    from hailtop.batch import ResourceGroup
    from hailtop.batch.job import Job
    from cpg_workflows.jobs.trim import trim
    from cpg_workflows.jobs.align_rna import align
    from cpg_workflows.jobs.markdups import markdup
    from cpg_workflows.jobs.bam_to_cram import bam_to_cram, cram_to_bam
    from cpg_workflows.jobs.count import count
    from cpg_workflows.jobs.outrider import outrider
    from cpg_workflows.jobs.fraser import fraser
    from cpg_workflows.filetypes import FastqPairs, FastqPair, BamPath, CramPath

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

    def capture_trim_cmd(*args, **kwargs) -> tuple[Job | None, FastqPair]:
        trim_job_output = trim(*args, **kwargs)
        j = trim_job_output[0]
        if j and isinstance(j, Job):
            cmd_str_list.append(
                '===== FASTQ TRIM JOB START =====\n\n' +
                '\n'.join(j._command) +
                '\n\n===== FASTQ TRIM JOB END =====\n\n'
            )
        return trim_job_output
    
    def capture_align_cmd(*args, **kwargs) -> tuple[list[Job], ResourceGroup] | tuple[None, BamPath]:
        align_job_output = align(*args, **kwargs)
        align_jobs = align_job_output[0]
        if align_jobs and isinstance(align_jobs, list) and all([isinstance(j, Job) for j in align_jobs]):
            cmd_str_list.append(
                '===== ALIGN STAGE START =====\n\n' +
                '----- Align sub-job start -----\n\n' +
                '\n\n----- Align sub-job end\n\n-----Align sub-job start -----\n\n'.join(['\n'.join(j._command) for j in align_jobs]) +
                '\n\n----- Align sub-job end -----\n\n'
                '\n\n===== ALIGN STAGE END =====\n\n'
            )
        return align_job_output
    
    def capture_markdup_cmd(*args, **kwargs) -> tuple[Job, ResourceGroup] | tuple[None, BamPath]:
        markdup_job_output = markdup(*args, **kwargs)
        markdup_job = markdup_job_output[0]
        if markdup_job:
            cmd_str_list.append(
                '===== MARKDUP JOB START =====\n\n' +
                '\n'.join(markdup_job._command) +
                '\n\n===== MARKDUP JOB END =====\n\n'
            )
        return markdup_job_output
    
    def capture_bam_to_cram_cmd(*args, **kwargs) -> tuple[Job, ResourceGroup] | tuple[None, CramPath]:
        bam_to_cram_job_output = bam_to_cram(*args, **kwargs)
        bam_to_cram_job = bam_to_cram_job_output[0]
        if bam_to_cram_job:
            cmd_str_list.append(
                '===== BAM TO CRAM JOB START =====\n\n' +
                '\n'.join(bam_to_cram_job._command) +
                '\n\n===== BAM TO CRAM JOB END =====\n\n'
            )
        return bam_to_cram_job_output
    
    def capture_cram_to_bam_cmd(*args, **kwargs) -> Job:
        cram_to_bam_job_output = cram_to_bam(*args, **kwargs)
        cram_to_bam_job = cram_to_bam_job_output[0]
        if cram_to_bam_job:
            cmd_str_list.append(
                '===== CRAM TO BAM JOB START =====\n\n' +
                '\n'.join(cram_to_bam_job._command) +
                '\n\n===== CRAM TO BAM JOB END =====\n\n'
            )
        return cram_to_bam_job_output
    
    def capture_count_cmd(*args, **kwargs) -> Job:
        count_job = count(*args, **kwargs)
        cmd_str_list.append(
            '===== COUNT STAGE START =====\n\n' +
            '\n'.join(count_job._command) +
            '\n\n===== COUNT STAGE END =====\n\n'
        )
        return count_job
    
    def capture_outrider_cmd(*args, **kwargs) -> Job:
        outrider_job = outrider(*args, **kwargs)
        cmd_str_list.append(
            '===== OUTRIDER STAGE START =====\n\n' +
            '\n'.join(outrider_job._command) +
            '\n\n===== OUTRIDER STAGE END =====\n\n'
        )
        return outrider_job
    
    def capture_fraser_cmd(*args, **kwargs) -> Job:
        fraser_job = fraser(*args, **kwargs)
        cmd_str_list.append(
            '===== FRASER STAGE START =====\n\n' +
            '\n'.join(fraser_job._command) +
            '\n\n===== FRASER STAGE END =====\n\n'
        )
        return fraser_job

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

    # Patch the trim function to capture the job command
    mocker.patch('cpg_workflows.jobs.trim.trim', capture_trim_cmd)
    # Patch the align function to capture the job command
    mocker.patch('cpg_workflows.jobs.align_rna.align', capture_align_cmd)
    # Patch the markdup function to capture the job command
    mocker.patch('cpg_workflows.jobs.markdups.markdup', capture_markdup_cmd)
    # Patch the bam_to_cram function to capture the job command
    mocker.patch('cpg_workflows.jobs.bam_to_cram.bam_to_cram', capture_bam_to_cram_cmd)
    # Patch the cram_to_bam function to capture the job command
    mocker.patch('cpg_workflows.jobs.bam_to_cram.cram_to_bam', capture_cram_to_bam_cmd)
    # Patch the count function to capture the job command
    mocker.patch('cpg_workflows.jobs.count.count', capture_count_cmd)
    # Patch the outrider function to capture the job calls
    mocker.patch('cpg_workflows.jobs.outrider.outrider', capture_outrider_cmd)
    # Patch the fraser function to capture the job calls
    mocker.patch('cpg_workflows.jobs.fraser.fraser', capture_fraser_cmd)

    from cpg_workflows.batch import get_batch
    from cpg_workflows.inputs import get_cohort
    from cpg_workflows.workflow import get_workflow

    # Imports specific for testing the trim stage
    from cpg_utils.hail_batch import dataset_path
    from cpg_workflows.stages.outrider import Outrider
    from cpg_workflows.stages.fraser import Fraser

    get_workflow().run(stages=[Outrider, Fraser])

    b = get_batch()
    trim_job = b.job_by_tool['fastp']
    align_job = b.job_by_tool['STAR']
    samtools_job = b.job_by_tool['samtools']
    markdup_job = b.job_by_tool['sambamba']
    bam_to_cram_job = b.job_by_tool['samtools_view']
    cram_to_bam_job = b.job_by_tool['samtools_view_cram_to_bam']
    featureCounts_job = b.job_by_tool['featureCounts']
    outrider_job = b.job_by_tool['outrider']
    fraser_job = b.job_by_tool['fraser']
    sample_list = get_cohort().get_sequencing_groups()

    # The number of FASTQ trim jobs should equal the number of FASTQ pairs
    alignment_input_list = [
        s.alignment_input_by_seq_type.get('transcriptome')
        for s in sample_list
    ]
    n_trim_jobs_list = [
        len(i) if isinstance(i, FastqPairs) else 0
        for i in alignment_input_list
    ]
    n_trim_jobs = sum(n_trim_jobs_list)
    assert trim_job['job_n'] == n_trim_jobs

    # The number of align jobs for each sample should equal:
    # the number of FASTQ pairs for that sample
    # + 1 if there are multiple FASTQ pairs (for merging the BAMs)
    # + 1 for the final sort and index job
    n_align_jobs_list = [
        len(i) + int(len(i) > 1) + 1 if isinstance(i, FastqPairs) else 0
        for i in alignment_input_list
    ]
    n_align_jobs = sum(n_align_jobs_list)
    assert align_job['job_n'] + samtools_job['job_n'] == n_align_jobs

    # The number of markdup jobs should equal the number of samples
    n_markdup_jobs = len(sample_list)
    assert markdup_job['job_n'] == n_markdup_jobs
    
    # The number of count jobs should equal the number of samples
    n_count_jobs = len(sample_list)
    assert featureCounts_job['job_n'] == n_count_jobs
    # The number of bam_to_cram jobs should equal the number of samples
    n_bam_to_cram_jobs = len(sample_list)
    assert bam_to_cram_job['job_n'] == n_bam_to_cram_jobs

    # The number of cram_to_bam jobs should equal the number of samples x2
    # for OUTRIDER and FRASER
    n_cram_to_bam_jobs = len(sample_list) * 2
    assert cram_to_bam_job['job_n'] == n_cram_to_bam_jobs

    # The number of outrider jobs should be 1 (for the cohort)
    n_outrider_jobs = 1
    assert outrider_job['job_n'] == n_outrider_jobs

    # The number of fraser jobs should be 1 (for the cohort)
    n_fraser_jobs = 1
    assert fraser_job['job_n'] == n_fraser_jobs

    output_path = dataset_path('cmd.txt')
    output_path_parent = os.path.dirname(output_path)
    makedirs(output_path_parent, exist_ok=True)
    with open(output_path, 'w') as f:
        for cmd_str in cmd_str_list:
            f.write(cmd_str)
        f.write('\n')
