"""
Test Hail Query functions.
"""
import shutil

import hail as hl
import pytest
import toml
from cpg_utils import to_path, Path
from cpg_utils.config import get_config, set_config_paths, update_dict
from cpg_utils.flows.batch import get_batch
from cpg_utils.flows.utils import timestamp
from cpg_utils.hail_batch import dataset_path, init_batch


DEFAULT_CONF = """
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
check_inputs = false
check_intermediates = false
check_expected_outputs = false
jc_intervals_num = 2
vqsr_intervals_num = 2
vep_intervals_num = 2
sequencing_type = 'genome'

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = true
"""


@pytest.fixture()
def tmp_dir() -> Path:
    dir_path = to_path('results') / timestamp()
    dir_path.mkdir(parents=True, exist_ok=True)
    yield dir_path
    shutil.rmtree(dir_path)


def _local_backend_config(dir_path: Path):
    _set_config(dir_path, {
        'workflow': {
            'path_scheme': 'local', 
            'local_dir': str(dir_path),
        },
        'hail': {'backend': 'local'},
    })


def _set_config(dir_path: Path, extra_conf: dict | None = None):
    d = toml.loads(DEFAULT_CONF)
    if extra_conf:
        update_dict(d, extra_conf)
    config_path = dir_path / 'config.toml'
    with config_path.open('w') as f:
        toml.dump(d, f)
    set_config_paths([str(config_path)])


def test_workflow_dry(tmp_dir: Path):
    """
    Run entire seqr-loader in a dry mode.
    """
    intervals_num = 2
    sample_ids = ['CPG56564', 'CPG56572']
    _set_config(tmp_dir, extra_conf={
        'workflow': {
            'stages': 'MtToEs',
            'only_samples': sample_ids,
            'skip_stages': ['Align'],
            'hc_intervals_num': intervals_num,
            'jc_intervals_num': intervals_num,
        },
        'hail': {
            'dry_run': True,
        }
    })
    import main
    main.main()
    assert get_batch().job_by_tool['gatk_HaplotypeCaller']['job_n'] == len(sample_ids) * intervals_num
    assert get_batch().job_by_tool['picard_MergeVcfs']['job_n'] == len(sample_ids)
    assert get_batch().job_by_tool['gatk_ReblockGVCF']['job_n'] == len(sample_ids)
    assert get_batch().job_by_tool['picard_CollectVariantCallingMetrics']['job_n'] == len(sample_ids) + 1
    assert get_batch().job_by_tool['gatk_GenomicsDBImport']['job_n'] == intervals_num


def test_annotate_cohort():
    """
    Test a query job.
    """
    chrom = 'chr20'
    locus1 = '5111495'
    locus2 = '5111607'
    interval = f'{chrom}-{locus1}-{locus2}'

    vcf_path = to_path(
        f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
        f'joint-called-{interval}.vcf.gz'
    )
    siteonly_vqsr_vcf_path = to_path(
        f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
        f'siteonly-vqsr-{interval}.vcf.gz'
    )
    vep_ht_path = to_path(
        f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/vep/{interval}.ht'
    )
    out_prefix = to_path(f'gs://cpg-fewgenomes-test/unittest/outputs')
    out_mt_path = out_prefix / f'cohort-{interval}.mt'

    from jobs.seqr_loader import annotate_cohort_jobs
    j = annotate_cohort_jobs(
        get_batch('test_annotate_cohort'),
        vcf_path=vcf_path,
        siteonly_vqsr_vcf_path=siteonly_vqsr_vcf_path,
        vep_ht_path=vep_ht_path,
        out_mt_path=out_mt_path,
        sequencing_type='genome',
        checkpoint_prefix=out_prefix / 'checkpoints',
    )
    get_batch().run(wait=True)

    # Testing
    init_batch()
    mt = hl.read_matrix_table(str(out_mt_path))
    mt.rows().show()
    assert mt.topmed.AC.collect() == [20555, 359, 20187]
    assert set(mt.geneIds.collect()[0]) == {'ENSG00000089063'}


def test_workflow():
    """
    Run entire workflow with a local backend.
    """
    # Local file system doesn't create directories automatically,
    # so creating the parent dataset path explicitly for all tests,
    # sort of simulating the cloud file system behaviour.
    to_path(dataset_path('')).parent.mkdir(parents=True, exist_ok=True)
    assert get_config()['workflow']['path_scheme'] == 'local'
