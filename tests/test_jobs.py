"""
Test individual jobs on Hail Batch.
"""

import shutil
import tempfile
import time
import unittest
from os.path import join
from typing import Dict

from hailtop.batch.job import Job
from analysis_runner import dataproc

from cpg_pipes import Path, to_path, Namespace
from cpg_pipes import benchmark, utils
from cpg_pipes.jobs import vep
from cpg_pipes.jobs.align import align
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.seqr_loader import annotate_dataset
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.types import CramPath, SequencingType
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes import images

try:
    from .utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID
except ImportError:
    from utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID  # type: ignore


def _read_file(path: Path) -> str:
    with path.open() as f:
        return f.read().strip()


class TestJobs(unittest.TestCase):
    """
    Test individual jobs
    """
    @property
    def out_bucket(self):
        return BASE_BUCKET / self.name / self.timestamp

    @property
    def tmp_bucket(self):
        return BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def setUp(self):
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.local_tmp_dir = tempfile.mkdtemp()

        self.pipeline = create_pipeline(
            name=self.name,
            description=self.name,
            analysis_dataset=DATASET,
            namespace=Namespace.TEST,
        )
        self.sequencing_type = SequencingType.WGS
        self.dataset = self.pipeline.add_dataset(DATASET)
        sample_name = f'Test-{self.timestamp}'
        self.sample = self.dataset.add_sample(sample_name, sample_name)
        self.sequencing_type = self.sequencing_type
        self.refs = self.pipeline.refs
        
    def tearDown(self) -> None:
        shutil.rmtree(self.local_tmp_dir)

    def _job_get_cram_details(
        self, 
        cram_path: Path, 
        out_bucket: Path,
        jobs: list[Job],
    ) -> Dict[str, Path]:
        """ 
        Add job that gets details of a CRAM file
        """
        test_j = self.pipeline.b.new_job('Parse CRAM sample name')
        test_j.image(images.SAMTOOLS_PICARD_IMAGE)
        sed = r's/.*SM:\([^\t]*\).*/\1/g'
        fasta_reference = self.refs.fasta_res_group(self.pipeline.b)
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram_path} '
            f'-c > {test_j.reads_num}'
        )
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram_path} '
            f'-c -f2 > {test_j.reads_num_mapped_in_proper_pair}'
        )
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram_path} '
            f'-H | grep \'^@RG\' | sed "{sed}" | uniq > {test_j.sample_name}'
        )
        test_j.depends_on(*jobs)
        d = {}
        for key in ['reads_num', 'reads_num_mapped_in_proper_pair', 'sample_name']:
            out_path = out_bucket / f'{key}.out'
            self.pipeline.b.write_output(test_j[key], str(out_path))
            d[key] = out_path
        return d

    def _job_get_gvcf_header(
        self, 
        vcf_path: Path, 
        jobs: list[Job],
    ) -> Path:
        """ 
        Parses header of GVCF file
        """
        vcf_input = self.pipeline.b.read_input(str(vcf_path))
        test_j = self.pipeline.b.new_job('Parse GVCF sample name')
        test_j.image(images.BCFTOOLS_IMAGE)
        test_j.command(f'bcftools query -l {vcf_input} > {test_j.output}')
        test_j.depends_on(*jobs)

        out_path = self.out_bucket / f'{self.sample.id}.out'
        self.pipeline.b.write_output(test_j.output, str(out_path))
        return out_path

    def test_alignment_fastq(self):
        """
        Test alignment job on a set of two FASTQ pairs 
        (tests processing in parallel merging and merging)
        """
        output_path = self.out_bucket / 'align_fastq' / 'result.cram'
        qc_bucket = self.out_bucket / 'align_fastq' / 'qc'
        
        jobs = align(
            b=self.pipeline.b,
            alignment_input=benchmark.tiny_fq,
            output_path=output_path,
            qc_bucket=qc_bucket,
            sample_name=self.sample.id,
            refs=self.refs,
        )
        cram_details_paths = self._job_get_cram_details(
            output_path, 
            out_bucket=self.out_bucket / 'align_fastq',
            jobs=jobs,
        )
        self.pipeline.submit_batch(wait=True)     
        
        self.assertTrue((
            qc_bucket / 
            'duplicate-metrics' / 
            f'{self.sample.id}-duplicate-metrics.csv'
        ).exists())

        self.assertEqual(
            self.sample.id, 
            _read_file(cram_details_paths['sample_name'])
        )
        self.assertAlmostEqual(
            20296, 
            int(_read_file(cram_details_paths['reads_num'])),
            delta=10,
        )
        self.assertAlmostEqual(
            18797,
            int(_read_file(cram_details_paths['reads_num_mapped_in_proper_pair'])),
            delta=10,
        )

    def test_alignment_cram(self):
        """
        Test alignment job on a CRAM input (tests realignment with bazam)
        """
        output_path = self.out_bucket / 'align_fastq' / 'result.cram'
        jobs = align(
            self.pipeline.b,
            alignment_input=benchmark.tiny_cram,
            output_path=output_path,
            sample_name=self.sample.id,
            refs=self.refs,
        )
        cram_details_paths = self._job_get_cram_details(
            output_path, 
            out_bucket=self.out_bucket / 'align',
            jobs=jobs,
        )
        self.pipeline.submit_batch(wait=True)

        self.assertEqual(
            self.sample.id, 
            _read_file(cram_details_paths['sample_name'])
        )
        self.assertAlmostEqual(
            285438,
            int(_read_file(cram_details_paths['reads_num'])),
            delta=50
        )
        self.assertAlmostEqual(
            273658,
            int(_read_file(cram_details_paths['reads_num_mapped_in_proper_pair'])),
            delta=10
        )

    def test_haplotype_calling(self):
        """
        Test individual sample haplotype calling
        """
        cram_path = CramPath(
            to_path(benchmark.BENCHMARK_BUCKET) /
            'outputs' / 'TINY_CRAM' / 'dragmap-picard.cram'
        )
        out_gvcf_path = self.out_bucket / 'test_haplotype_calling.g.vcf.gz'
        jobs = produce_gvcf(
            self.pipeline.b,
            sample_name=self.sample.id,
            refs=self.refs,
            cram_path=cram_path,
            number_of_intervals=2,
            tmp_bucket=self.tmp_bucket,
            dragen_mode=True,
            output_path=out_gvcf_path,
            sequencing_type=self.sequencing_type,
        )
        test_result_path = self._job_get_gvcf_header(out_gvcf_path, jobs)
        self.pipeline.submit_batch(wait=True)     
        contents = _read_file(test_result_path)
        self.assertEqual(self.sample.id, contents.split()[-1])

    def test_joint_calling(self):
        """
        Test joint variant calling
        """
        genomicsdb_bucket = self.out_bucket / 'genomicsdb'
        out_vcf_path = self.out_bucket / 'joint-called.vcf.gz'
        out_siteonly_vcf_path = to_path(
            str(out_vcf_path).replace('.vcf.gz', '-siteonly.vcf.gz')
        )

        proj = self.pipeline.cohort.add_dataset(DATASET)
        for sid in SAMPLES:
            proj.add_sample(sid, sid)

        jobs = make_joint_genotyping_jobs(
            b=self.pipeline.b,
            out_vcf_path=out_vcf_path,
            out_siteonly_vcf_path=out_siteonly_vcf_path,
            refs=self.refs,
            samples=self.pipeline.get_all_samples(),
            genomicsdb_bucket=genomicsdb_bucket,
            gvcf_by_sid=SUBSET_GVCF_BY_SID,
            tmp_bucket=self.tmp_bucket,
            overwrite=True,
            scatter_count=10,
            tool=JointGenotyperTool.GenotypeGVCFs,
            sequencing_type=self.pipeline.cohort.get_sequencing_type(),
        )
        test_result_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.exists(out_vcf_path))
        self.assertTrue(utils.exists(out_siteonly_vcf_path))
        contents = _read_file(test_result_path)
        self.assertEqual(len(SAMPLES), len(contents.split()))
        self.assertEqual(set(SAMPLES), set(contents.split()))

    def test_vqsr(self):
        """
        Test AS-VQSR
        """
        siteonly_vcf_path = to_path(
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'joint-called-siteonly.vcf.gz'
        )
        tmp_vqsr_bucket = self.tmp_bucket / 'vqsr'
        out_vcf_path = self.out_bucket / 'vqsr' / 'vqsr.vcf.gz'
        jobs = make_vqsr_jobs(
            b=self.pipeline.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            refs=self.refs,
            work_bucket=tmp_vqsr_bucket,
            gvcf_count=len(SAMPLES),
            scatter_count=10,
            output_vcf_path=out_vcf_path,
            use_as_annotations=True,
            overwrite=True,
            sequencing_type=self.sequencing_type,
        )
        res_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.exists(out_vcf_path))
        contents = _read_file(res_path)
        self.assertEqual(0, len(contents.split()))  # site-only doesn't have any samples

    def test_vep(self):
        """
        Tests command line VEP with LoF plugin and MANE_SELECT annotation
        """
        self.timestamp = '20220326-1639'

        site_only_vcf_path = to_path(
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'vqsr.vcf.gz'
        )
        out_vcf_path = self.out_bucket / 'vep' / 'vep.vcf.gz'
        jobs = vep.vep(
            self.pipeline.b, 
            vcf_path=site_only_vcf_path,
            refs=self.refs,
            sequencing_type=self.sequencing_type,
            out_vcf_path=out_vcf_path,
            scatter_count=10,
            overwrite=False,
        )
        
        # Add test job
        test_j = self.pipeline.b.new_job('Parse GVCF sample name')
        test_j.image(images.BCFTOOLS_IMAGE)
        test_j.command(f"""
        bcftools +split-vep {self.pipeline.b.read_input(str(out_vcf_path))} \
        -f '%CHROM:%POS %SYMBOL %BIOTYPE %MANE_SELECT %LoF %LoF_filter\n' \
        -i'BIOTYPE="protein_coding" & LoF_filter="ANC_ALLELE"' -s worst \
        > {test_j.output}
        """)
        test_j.depends_on(*jobs)
        test_out_path = self.out_bucket / f'{self.sample.id}.out'
        self.pipeline.b.write_output(test_j.output, str(test_out_path))
        
        # Run Batch
        self.pipeline.submit_batch(wait=True)
        
        # Check results
        self.assertTrue(out_vcf_path.exists())
        contents = _read_file(test_out_path)
        self.assertEqual(
            'chr20:5111495 TMEM230 protein_coding NM_001009923.2 LC ANC_ALLELE', 
            contents
        )

    def test_seqr_loader_annotate_dataset_seqr(self):
        """
        Test subset_mt.py script from seqr_loader.
        """
        dataset = 'seqr'
        pipeline = create_pipeline(
            name=self.name,
            description=self.name,
            analysis_dataset=dataset,
            namespace=Namespace.MAIN,
        )
        cohort_mt_path = to_path('gs://cpg-seqr-main-tmp/mt/combined.mt')
        dataset_mt_path = to_path(f'gs://cpg-circa-main-analysis/mt/circa-{self.timestamp}.mt')
        samples = ['CPG200675', 'CPG200642', 'CPG200501', 'CPG200527', 'CPG200519', 'CPG200410', 'CPG200477', 'CPG200485', 'CPG200493', 'CPG200444', 'CPG200428', 'CPG200451', 'CPG200436', 'CPG200659', 'CPG200667', 'CPG200535', 'CPG200543', 'CPG200550', 'CPG200584', 'CPG200568', 'CPG200592', 'CPG200600', 'CPG200626', 'CPG200618', 'CPG200683']
        tmp_bucket = to_path('gs://cpg-seqr-main-tmp')
        annotate_dataset(
            b=pipeline.b,
            annotated_mt_path=cohort_mt_path,
            sample_ids=samples,
            tmp_bucket=tmp_bucket,
            output_mt_path=dataset_mt_path,
            hail_billing_project=dataset,
            hail_bucket=tmp_bucket,
        )
        pipeline.submit_batch(wait=False)     

    def test_seqr_loader_annotate_dataset(self):
        """
        Test subset_mt.py script from seqr_loader.
        """
        cohort_mt_path = to_path(
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/seqr_loader/cohort.mt'
        )
        dataset_mt_path = self.out_bucket / 'dataset.mt'
        annotate_dataset(
            b=self.pipeline.b,
            annotated_mt_path=cohort_mt_path,
            sample_ids=SAMPLES[:3],
            tmp_bucket=self.tmp_bucket,
            output_mt_path=dataset_mt_path,
            hail_billing_project=DATASET,
            hail_bucket=self.tmp_bucket,
        )
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.exists(dataset_mt_path))

    def test_seqr_loader(self):
        """
        Assuming variants are called and VQSR'ed, tests loading
        into a matrix table and annotation for Seqr
        """
        vcf_path = (
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'joint-called.vcf.gz'
        )
        vqsr_vcf_path = (
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'vqsr.vcf.gz'
        )
        
        cohort_mt_path = f'{self.out_bucket}/seqr_loader/cohort.mt'
        dataset_mt_path = f'{self.out_bucket}/seqr_loader/dataset.mt'
        dataset_vcf_path = f'{self.out_bucket}/seqr_loader/dataset.vcf.bgz'

        cluster = dataproc.setup_dataproc(
            self.pipeline.b,
            cluster_name='Test seqr loader',
            max_age='1h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            # default VEP initialization script (used with --vep) installs VEP=v95,
            # but we need v105, so we use a modified vep-GRCh38.sh (with --init) from
            # this repo: production-pipelines/vep/vep-GRCh38.sh
            init=['gs://cpg-reference/vep/vep-GRCh38.sh'], 
            scopes=['cloud-platform'],
        )
        vcf_to_mt_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "vcf_to_mt.py")} '
            f'--vcf-path {vcf_path} '
            f'--site-only-vqsr-vcf-path {vqsr_vcf_path} '
            f'--dest-mt-path {cohort_mt_path} '
            f'--bucket {self.out_bucket}/seqr_loader '
            f'--disable-validation '
            f'--make-checkpoints',
            job_name='Make MT and annotate cohort',
        )
        mt_to_datasetmt_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "subset_mt.py")} '
            f'--mt-path {cohort_mt_path} '
            f'--out-mt-path {dataset_mt_path}',
            job_name=f'Annotate dataset',
        )
        mt_to_datasetmt_j.depends_on(vcf_to_mt_j)
        datasetmt_to_es_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_es.py")} '
            f'--mt-path {dataset_mt_path} '
            f'--es-index test-{self.timestamp} '
            f'--es-index-min-num-shards 1',
            job_name=f'Create ES index',
        )
        datasetmt_to_es_j.depends_on(mt_to_datasetmt_j)
        datasetmt_to_vcf_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_vcf.py")} '
            f'--mt-path {dataset_mt_path} '
            f'--out-vcf-path {dataset_vcf_path}',
            job_name=f'Convert to VCF',
        )
        datasetmt_to_vcf_j.depends_on(mt_to_datasetmt_j)
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.exists(dataset_vcf_path))
