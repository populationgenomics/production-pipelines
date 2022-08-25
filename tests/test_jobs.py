"""
Test individual jobs on Hail Batch.
"""
import logging
import shutil
import tempfile
import unittest
from typing import Dict
from unittest import skip

from cpg_utils.flows.utils import timestamp
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group, init_batch
from cpg_utils.flows.batch import setup_batch
from cpg_utils.flows.targets import Dataset
from cpg_utils.flows.filetypes import FastqPair, CramPath, FastqPairs

from jobs import vep
from jobs.align import align
from jobs.haplotype_caller import produce_gvcf
from jobs.joint_genotyping import (
    make_joint_genotyping_jobs,
    JointGenotyperTool,
)
from jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from jobs.somalier import check_pedigree_job
from jobs.vqsr import make_vqsr_jobs

try:
    import utils
except ModuleNotFoundError:
    from . import utils


logger = logging.getLogger(__file__)


BENCHMARK_BUCKET = to_path('gs://cpg-fewgenomes-test/benchmark')
TOY_INPUTS_BUCKET = BENCHMARK_BUCKET / 'inputs/toy'
RESULTS_BUCKET = BENCHMARK_BUCKET / 'outputs'


# 40k reads:
tiny_fq = FastqPairs(
    [
        # This set is 50MB each:
        FastqPair(
            TOY_INPUTS_BUCKET / '2-699835.L001.R1.n40000.fastq.gz',
            TOY_INPUTS_BUCKET / '2-699835.L001.R2.n40000.fastq.gz',
        ),
        FastqPair(
            TOY_INPUTS_BUCKET / '2-699835.L002.R1.n40000.fastq.gz',
            TOY_INPUTS_BUCKET / '2-699835.L002.R2.n40000.fastq.gz',
        ),
    ]
)

# ~300k reads:
tiny_cram = CramPath(TOY_INPUTS_BUCKET / 'NA12878-chr21-tiny.cram')

# WGS:
giab_crams = {
    sn: CramPath(f'gs://cpg-reference/validation/giab/cram/{sn}.cram')
    for sn in ['NA12878', 'NA12891', 'NA12892']
}
na12878fq = FastqPairs(
    [
        FastqPair(
            BENCHMARK_BUCKET / 'inputs/NA12878/ERR194147_1.fastq.gz',
            BENCHMARK_BUCKET / 'inputs/NA12878/ERR194147_2.fastq.gz',
        )
    ]
)
perth_neuro_fq = FastqPairs(
    [
        FastqPair(
            BENCHMARK_BUCKET
            / 'inputs/PERTHNEURO_FQ/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R1.fastq.gz',
            BENCHMARK_BUCKET
            / 'inputs/PERTHNEURO_FQ/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R2.fastq.gz',
        )
    ]
)
perth_neuro_cram = CramPath(
    BENCHMARK_BUCKET / 'inputs/PERTHNEURO_CRAM/CPG13045.cram',
)


def _read_file(path: Path) -> str:
    with path.open() as f:
        return f.read().strip()


class TestJobs(unittest.TestCase):
    """
    Test individual jobs
    """

    @property
    def out_bucket(self):
        """Write test results here"""
        return utils.BASE_BUCKET / self.name / self.timestamp

    @property
    def tmp_bucket(self):
        """Use for intermediate/tmp outputs"""
        return utils.BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def setUp(self):
        """Set up test environment"""
        self.name = self._testMethodName
        self.timestamp = timestamp()
        logger.info(f'Timestamp: {self.timestamp}')
        self.local_tmp_dir = tempfile.mkdtemp()
        utils.setup_env(self.tmp_bucket)

        self.batch = setup_batch(self.name)
        self.dataset = Dataset(utils.DATASET)
        sample_name = f'Test-{self.timestamp}'
        self.sample = self.dataset.add_sample(sample_name, sample_name)

        # Interval to take on chr20:
        self.chrom = 'chr20'
        self.locus1 = '5111495'
        self.locus2 = '5111607'
        self.interval = f'{self.chrom}-{self.locus1}-{self.locus2}'

    def tearDown(self) -> None:
        """Remove test environment"""
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
        test_j = self.batch.new_job('Parse CRAM sample name')
        test_j.image(image_path('samtools'))
        sed = r's/.*SM:\([^\t]*\).*/\1/g'
        fasta_reference = fasta_res_group(self.batch)
        cram = self.batch.read_input_group(
            **{
                'cram': str(cram_path),
                'cram.crai': str(cram_path) + '.crai',
            }
        )
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram.cram} '
            f'-c > {test_j.reads_num}'
        )
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram.cram} '
            f'-c -f2 > {test_j.reads_num_mapped_in_proper_pair}'
        )
        test_j.command(
            f'samtools view -T {fasta_reference.base} {cram.cram} '
            f'-H | grep \'^@RG\' | sed "{sed}" | uniq > {test_j.sample_name}'
        )
        test_j.depends_on(*jobs)
        d = {}
        for key in ['reads_num', 'reads_num_mapped_in_proper_pair', 'sample_name']:
            out_path = out_bucket / f'{key}.out'
            self.batch.write_output(test_j[key], str(out_path))
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
        vcf_input = self.batch.read_input(str(vcf_path))
        test_j = self.batch.new_job('Parse GVCF sample name')
        test_j.image(image_path('bcftools'))
        test_j.command(f'bcftools query -l {vcf_input} > {test_j.output}')
        test_j.depends_on(*jobs)

        out_path = self.out_bucket / f'{self.sample.id}.out'
        self.batch.write_output(test_j.output, str(out_path))
        return out_path

    def test_align_fastq(self):
        """
        Test alignment job on a set of two FASTQ pairs
        (tests processing in parallel merging and merging)
        """
        output_path = self.out_bucket / 'align_fastq' / 'result.cram'
        self.sample.alignment_input_by_seq_type['genome'] = tiny_fq
        jobs = align(
            b=self.batch,
            sample=self.sample,
            output_path=output_path,
        )
        cram_details_paths = self._job_get_cram_details(
            output_path,
            out_bucket=self.out_bucket / 'align_fastq',
            jobs=jobs,
        )
        self.batch.run(wait=False)
        self.assertEqual(self.sample.id, _read_file(cram_details_paths['sample_name']))
        self.assertAlmostEqual(
            20180,
            int(_read_file(cram_details_paths['reads_num'])),
            delta=10,
        )
        self.assertAlmostEqual(
            14536,
            int(_read_file(cram_details_paths['reads_num_mapped_in_proper_pair'])),
            delta=10,
        )

    def test_align_cram(self):
        """
        Test alignment job on a CRAM input (tests realignment with bazam)
        """
        sid = 'CPG196519'
        cram_path = CramPath(utils.BASE_BUCKET / f'inputs/toy/cram/{sid}.cram')
        self.sample.alignment_input_by_seq_type['genome'] = cram_path
        output_path = self.out_bucket / 'result.cram'
        jobs = align(
            self.batch,
            sample=self.sample,
            output_path=output_path,
            realignment_shards_num=4,
        )
        cram_details_paths = self._job_get_cram_details(
            output_path,
            out_bucket=self.out_bucket / 'align',
            jobs=jobs,
        )
        self.batch.run(wait=False)

        self.assertEqual(sid, _read_file(cram_details_paths['sample_name']))
        cram_details = {k: _read_file(v) for k, v in cram_details_paths.items()}
        print(cram_details_paths)
        self.assertAlmostEqual(223007, int(cram_details['reads_num']), delta=100)
        self.assertAlmostEqual(
            222678,
            int(cram_details['reads_num_mapped_in_proper_pair']),
            delta=100,
        )

    def test_genotype(self):
        """
        Test individual sample haplotype calling.
        """
        sid = 'CPG196519'
        cram_path = CramPath(
            utils.BASE_BUCKET / f'inputs/toy/cram_realigned/{sid}.cram'
        )
        out_gvcf_path = self.out_bucket / 'test.g.vcf.gz'
        jobs = produce_gvcf(
            self.batch,
            sample_name=sid,
            cram_path=cram_path,
            scatter_count=4,
            tmp_prefix=self.tmp_bucket,
            dragen_mode=True,
            output_path=out_gvcf_path,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        test_result_path = self._job_get_gvcf_header(out_gvcf_path, jobs)
        self.batch.run(wait=True)
        contents = _read_file(test_result_path)
        self.assertEqual(sid, contents.split()[-1])

    def test_joint_calling(self):
        """
        Test joint variant calling.
        """
        out_vcf_path = self.out_bucket / 'joint-called.vcf.gz'
        out_siteonly_vcf_path = to_path(
            str(out_vcf_path).replace('.vcf.gz', '-siteonly.vcf.gz')
        )

        ds = Dataset(utils.DATASET)
        for sid in utils.SAMPLES:
            ds.add_sample(sid, sid)

        jobs = make_joint_genotyping_jobs(
            b=self.batch,
            out_vcf_path=out_vcf_path,
            out_siteonly_vcf_path=out_siteonly_vcf_path,
            gvcf_by_sid=utils.EXOME_1PCT_GVCF_BY_SID,
            tmp_bucket=self.tmp_bucket,
            overwrite=True,
            scatter_count=4,
            tool=JointGenotyperTool.GenotypeGVCFs,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        test_result_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.batch.run(wait=True)
        self.assertTrue(out_vcf_path.exists())
        self.assertTrue(out_siteonly_vcf_path.exists())
        contents = _read_file(test_result_path)
        self.assertEqual(len(utils.SAMPLES), len(contents.split()))
        self.assertEqual(set(utils.SAMPLES), set(contents.split()))

    def test_vqsr(self):
        """
        Test AS-VQSR. Needs 5% exome to avoid issues like
         "One or more annotations (usually MQ) may have insufficient variance.":
        https://batch.hail.populationgenomics.org.au/batches/55673/jobs/3
        """
        siteonly_vcf_path = (
            utils.BASE_BUCKET / 'inputs/exome5pct/9samples-joint-called-siteonly.vcf.gz'
        )
        tmp_vqsr_bucket = self.tmp_bucket / 'vqsr'
        out_vcf_path = self.out_bucket / 'vqsr' / 'vqsr.vcf.gz'
        jobs = make_vqsr_jobs(
            b=self.batch,
            input_vcf_or_mt_path=siteonly_vcf_path,
            tmp_prefix=tmp_vqsr_bucket,
            gvcf_count=len(utils.SAMPLES),
            scatter_count=4,
            out_path=out_vcf_path,
            use_as_annotations=True,
            overwrite=True,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome5pct/calling_regions.interval_list',
        )
        res_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.batch.run(wait=True)
        self.assertTrue(out_vcf_path.exists())
        contents = _read_file(res_path)
        self.assertEqual(0, len(contents.split()))  # site-only doesn't have any samples

    def test_vep_vcf(self):
        """
        Test VEP on VCF, parallelised by interval, with LOFtee plugin
        and MANE_SELECT annotation.
        """
        siteonly_vcf_path = (
            utils.BASE_BUCKET
            / 'inputs/exome0.1pct/9samples-joint-called-siteonly.vcf.gz'
        )
        out_vcf_path = self.out_bucket / 'vep' / 'vep.vcf.gz'
        jobs = vep.vep_jobs(
            self.batch,
            vcf_path=siteonly_vcf_path,
            out_path=out_vcf_path,
            scatter_count=4,
            overwrite=False,
            tmp_prefix=self.tmp_bucket,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )

        # Add test job
        test_j = self.batch.new_job('Parse VEP VCF')
        test_j.image(image_path('bcftools'))
        test_j.command(
            f"""
        bcftools +split-vep {self.batch.read_input(str(out_vcf_path))} \
        -f '%CHROM:%POS %SYMBOL %BIOTYPE %MANE_SELECT %LoF %LoF_filter\n' \
        -i'BIOTYPE="protein_coding" & LoF_filter="ANC_ALLELE"' -s worst \
        > {test_j.output}
        """
        )
        test_j.depends_on(*jobs)
        test_out_path = self.out_bucket / f'{self.sample.id}.out'
        self.batch.write_output(test_j.output, str(test_out_path))

        # Run Batch
        self.batch.run(wait=True)

        # Check results
        self.assertTrue(out_vcf_path.exists())
        contents = _read_file(test_out_path)
        self.assertListEqual(
            'chr8:22163875 SFTPC protein_coding . LC ANC_ALLELE'.split(),
            contents.split(),
        )

    def test_vep_to_json(self):
        """
        Tests a single VEP job that outputs JSON.
        """
        # Interval to take on chr20:
        locus1 = '5111495'
        locus2 = '5111607'
        siteonly_vqsr_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'siteonly-vqsr-chr20-{locus1}-{locus2}.vcf.gz'
        )
        out_path = self.out_bucket / 'vep_json' / f'chr20-{locus1}-{locus2}.json_list'
        j = vep.vep_one(
            self.batch,
            vcf=siteonly_vqsr_path,
            out_path=out_path,
            out_format='json',
        )

        # Add test job
        mane_transcript = 'NM_001009923.2'  # expected transcript on locus1
        test_j = self.batch.new_job('Parse VEP results')
        test_j.image(image_path('hail'))
        test_j.command(
            f"""
        cat {self.batch.read_input(str(out_path))} | zgrep {locus1} | \
        jq -r '.transcript_consequences[] | select(.mane_select=="{mane_transcript}") | .mane_select' \
        > {test_j.output}
        """
        )
        test_j.depends_on(j)
        test_out_path = self.out_bucket / f'{self.sample.id}.out'
        self.batch.write_output(test_j.output, str(test_out_path))

        # Run Batch
        self.batch.run(wait=True)

        # Check results
        self.assertTrue(out_path.exists())
        contents = _read_file(test_out_path)
        self.assertEqual(mane_transcript, contents)

    def test_vep_json_to_ht(self):
        """
        Test parsing VEP JSON into a Hail table.
        """
        vep_json_list_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/vep/'
            f'{self.interval}.json_list'
        )
        out_path = self.out_bucket / 'vep_ht' / f'{self.interval}.ht'

        vep.gather_vep_json_to_ht(
            b=self.batch,
            vep_results_paths=[vep_json_list_path],
            out_path=out_path,
        )
        self.batch.run(wait=True)

        import hail as hl

        init_batch()
        ht = hl.read_table(str(out_path))
        interval = hl.parse_locus_interval(
            f'{self.chrom}:{self.locus1}-{int(self.locus1) + 1}'
        )
        ht = hl.filter_intervals(ht, [interval])
        transcript = ht.vep.transcript_consequences.mane_select.collect()[0][1]
        lof = ht.vep.transcript_consequences.lof.collect()[0][1]
        self.assertEqual(lof, 'LC')
        self.assertEqual(transcript, 'NM_001009923.2')

    @skip('Large test')
    def test_seqr_loader_annotate_cohort_large(self):
        vcf_path = to_path(
            utils.BASE_BUCKET / 'inputs/exome5pct/9samples-joint-called.vcf.gz'
        )
        siteonly_vqsr_path = (
            utils.BASE_BUCKET
            / 'inputs/exome5pct/9samples-joint-called-siteonly-vqsr.vcf.gz'
        )
        vep_ht_path = to_path(utils.BASE_BUCKET / 'inputs/exome5pct/9samples-vep.ht')
        out_mt_path = self.out_bucket / f'cohort-exome5pct.mt'
        annotate_cohort_jobs(
            self.batch,
            vcf_path=vcf_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_path,
            vep_ht_path=vep_ht_path,
            output_mt_path=out_mt_path,
            checkpoint_prefix=self.tmp_bucket / 'checkpoints',
            sequencing_type='genome',
        )
        self.batch.run(wait=True)

    def test_seqr_loader_annotate_cohort(self):
        vcf_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'joint-called-{self.interval}.vcf.gz'
        )
        siteonly_vqsr_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'siteonly-vqsr-{self.interval}.vcf.gz'
        )
        vep_ht_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/' f'vep/{self.interval}.ht'
        )
        out_mt_path = self.out_bucket / f'cohort-{self.interval}.mt'
        annotate_cohort_jobs(
            self.batch,
            vcf_path=vcf_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_path,
            vep_ht_path=vep_ht_path,
            output_mt_path=out_mt_path,
            checkpoint_prefix=self.tmp_bucket / 'checkpoints',
            sequencing_type='genome',
        )
        self.batch.run(wait=True)

        # Testing
        import hail as hl

        init_batch()
        mt = hl.read_matrix_table(str(out_mt_path))
        mt.rows().show()
        self.assertListEqual(mt.topmed.AC.collect(), [20555, 359, 20187])
        self.assertSetEqual(set(mt.geneIds.collect()[0]), {'ENSG00000089063'})

    def test_seqr_loader_annotate_dataset(self):
        """
        Test subset_mt.py script from seqr_loader.
        """
        mt_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'seqr_loader/cohort-{self.interval}.mt'
        )
        out_mt_path = self.out_bucket / 'seqr_loader' / f'dataset-{self.interval}.mt'
        annotate_dataset_jobs(
            b=self.batch,
            mt_path=mt_path,
            sample_ids=utils.SAMPLES[:3],
            tmp_bucket=self.tmp_bucket,
            output_mt_path=out_mt_path,
        )
        self.batch.run(wait=True)

        # Testing
        self.assertTrue(out_mt_path.exists())
        import hail as hl

        init_batch()
        mt = hl.read_matrix_table(str(out_mt_path))
        mt.rows().show()
        self.assertListEqual(mt.s.collect(), utils.SAMPLES[:3])
        self.assertSetEqual(
            set(mt.samples_gq['20_to_25'].collect()[0]), {'CPG196519', 'CPG196527'}
        )
        self.assertSetEqual(set(mt.samples_ab['40_to_45'].collect()[0]), {'CPG196535'})

    def test_check_pedigree(self):
        inputs_bucket = utils.BASE_BUCKET / 'inputs' / 'check_pedigree'
        check_pedigree_job(
            self.batch,
            samples_file=self.batch.read_input(
                str(inputs_bucket / 'somalier-samples.tsv')
            ),
            pairs_file=self.batch.read_input(str(inputs_bucket / 'somalier-pairs.tsv')),
            expected_ped=self.batch.read_input(
                str(inputs_bucket / 'somalier-samples.tsv')
            ),
            dataset_name='fewgenomes',
        )
        self.batch.run(wait=True)
