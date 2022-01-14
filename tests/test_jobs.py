import shutil
import tempfile
import time
import unittest
from os.path import join, basename

# Make sure utils are imported beforehand, so the environemnt
# is set before cpg_pipes are imported
from .utils import BASE_BUCKET, PROJECT, SAMPLES, SUBSET_GVCF_BY_SID

from analysis_runner import dataproc

from cpg_pipes import benchmark
from cpg_pipes.jobs.align import align, Aligner
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.jobs.joint_genotyping import make_joint_genotyping_jobs, \
    JointGenotyperTool
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.pipeline import Pipeline, Sample
from cpg_pipes import utils, resources


class TestJobs(unittest.TestCase):
    """
    Test individual jobs
    """
    def setUp(self):
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.out_bucket = f'{BASE_BUCKET}/{self.name}/{self.timestamp}'
        self.tmp_bucket = f'{self.out_bucket}/tmp'
        self.local_tmp_dir = tempfile.mkdtemp()

        self.pipeline = Pipeline(
            name=self._testMethodName,
            title=self._testMethodName,
            analysis_project=PROJECT,
            output_version='v0',
            namespace='test',
        )
        self.sample_name = f'Test-{self.timestamp}'
        
    def tearDown(self) -> None:
        shutil.rmtree(self.local_tmp_dir)

    def _read_object_contents(self, path):
        local_tmp_dir = self.local_tmp_dir
        local_path = join(local_tmp_dir, basename(path))
        utils.gsutil_cp(path, local_path)
        with open(local_path) as f:
            return f.read().strip()

    def _job_get_cram_sample_name(self, cram_path) -> str:
        """ 
        Add job that gets sample name from CRAM file
        """
        test_j = self.pipeline.b.new_job('Parse CRAM sample name')
        test_j.image(resources.SAMTOOLS_PICARD_IMAGE)
        sed = r's/.*SM:\([^\t]*\).*/\1/g'
        test_j.command(
            f'samtools view -H {cram_path}'
            f' | grep \'^@RG\' | sed "{sed}" | uniq > {test_j.output}'
        )
        out_path = join(self.out_bucket, f'{self.sample_name}.out')
        self.pipeline.b.write_output(test_j.output, out_path)
        return out_path

    def _job_get_gvcf_header(self, gvcf_path) -> str:
        """ 
        Parses header of GVCF file
        """
        test_j = self.pipeline.b.new_job('Parse GVCF sample name')
        test_j.image(resources.BCFTOOLS_IMAGE)
        test_j.command(
            f'bcftools view {gvcf_path}'
            f' | awk \'/^#CHROM/ {{ print; exit }}\' > {test_j.output}'
        )
        out_path = join(self.out_bucket, f'{self.sample_name}.out')
        self.pipeline.b.write_output(test_j.output, out_path)
        return out_path

    @unittest.skip('Skip')
    def test_alignment(self):
        """
        Test alignment job
        TODO: use Stage objects. Add mocks to connect stages.
        """
        for inp in [benchmark.tiny_fq, benchmark.tiny_cram]:
            with self.subTest(inp=inp):
                j = align(
                    self.pipeline.b,
                    alignment_input=inp,
                    sample_name=self.sample_name,
                    project_name=PROJECT,
                    aligner=Aligner.DRAGMAP,
                )
                test_result_path = self._job_get_cram_sample_name(j.output_cram.cram)
                self.pipeline.submit_batch(wait=True)     
                contents = self._read_object_contents(test_result_path)
                self.assertEqual(self.sample_name, contents)

    @unittest.skip('Skip')
    def test_haplotype_calling(self):
        """
        Test individual sample haplotype calling
        """
        cram_path = f'{benchmark.BENCHMARK_BUCKET}/outputs/TINY_CRAM/dragmap-picard.cram'
        j = produce_gvcf(
            self.pipeline.b,
            sample_name=self.sample_name,
            project_name=PROJECT,
            cram_path=cram_path,
            number_of_intervals=10,
            tmp_bucket=self.tmp_bucket,
            dragen_mode=True,
        )
        test_result_path = self._job_get_gvcf_header(j.output_gvcf['g.vcf.gz'])
        self.pipeline.submit_batch(wait=True)     
        contents = self._read_object_contents(test_result_path)
        self.assertEqual(self.sample_name, contents.split()[-1])

    @unittest.skip('Skip')
    def test_joint_calling(self):
        """
        Test joint variant calling
        """
        genomicsdb_bucket = f'{self.out_bucket}/genomicsdb'
        out_vcf_path = f'{self.out_bucket}/joint-called.vcf.gz'
        out_siteonly_vcf_path = out_vcf_path.replace('.vcf.gz', '-siteonly.vcf.gz')

        j = make_joint_genotyping_jobs(
            b=self.pipeline.b,
            out_vcf_path=out_vcf_path,
            out_siteonly_vcf_path=out_siteonly_vcf_path,
            samples=[Sample(s, s) for s in SAMPLES],
            genomicsdb_bucket=genomicsdb_bucket,
            gvcf_by_sid=SUBSET_GVCF_BY_SID,
            tmp_bucket=self.tmp_bucket,
            overwrite=True,
            local_tmp_dir=self.local_tmp_dir,
            scatter_count=10,
            tool=JointGenotyperTool.GenotypeGVCFs,
        )
        test_result_path = self._job_get_gvcf_header(j.output_vcf['vcf.gz'])
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.file_exists(out_vcf_path))
        self.assertTrue(utils.file_exists(out_siteonly_vcf_path))
        contents = self._read_object_contents(test_result_path)
        self.assertEqual(9 + len(SAMPLES), len(contents.split()))

    @unittest.skip('Skip')
    def test_vqsr(self):
        """
        Test AS-VQSR
        """
        siteonly_vcf_path = (
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'joint-called-siteonly.vcf.gz'
        )
        tmp_vqsr_bucket = f'{self.tmp_bucket}/vqsr'
        out_vcf_path = f'{self.out_bucket}/vqsr/vqsr.vcf.gz'
        j = make_vqsr_jobs(
            b=self.pipeline.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            work_bucket=tmp_vqsr_bucket,
            gvcf_count=len(SAMPLES),
            scatter_count=10,
            output_vcf_path=out_vcf_path,
            use_as_annotations=True,
            overwrite=True,
        )
        res_path = self._job_get_gvcf_header(j.output_vcf['vcf.gz'])
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.file_exists(out_vcf_path))
        contents = self._read_object_contents(res_path)
        self.assertEqual(8, len(contents.split()))  # site-only doesn't have any samples
    
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
        project_mt_path = f'{self.out_bucket}/seqr_loader/project.mt'
        
        cluster = dataproc.setup_dataproc(
            self.pipeline.b,
            cluster_name='Test seqr loader',
            max_age='1h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            vep='GRCh38',            
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
        mt_to_projectmt_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "mt_to_projectmt.py")} '
            f'--mt-path {cohort_mt_path} '
            f'--out-mt-path {project_mt_path}',
            job_name=f'Annotate project',
        )
        mt_to_projectmt_j.depends_on(vcf_to_mt_j)
        projectmt_to_es_j = cluster.add_job(
            f'{join("..", utils.QUERY_SCRIPTS_DIR, "seqr", "projectmt_to_es.py")} '
            f'--mt-path {project_mt_path} '
            f'--es-index test-{self.timestamp} '
            f'--es-index-min-num-shards 1',
            job_name=f'Create ES index',
        )
        projectmt_to_es_j.depends_on(mt_to_projectmt_j)
        self.pipeline.submit_batch(wait=True)     
        self.assertTrue(utils.file_exists(project_mt_path))
