"""
Test individual jobs on Hail Batch.
"""
import logging
import shutil
import tempfile
import unittest
from typing import Dict
from unittest import skip

from hailtop.batch.job import Job

from cpg_pipes import Path, to_path, Namespace
from cpg_pipes import benchmark
from cpg_pipes import images
from cpg_pipes.hailquery import init_batch
from cpg_pipes.jobs import vep
from cpg_pipes.jobs.align import align
from cpg_pipes.jobs.haplotype_caller import produce_gvcf
from cpg_pipes.jobs.joint_genotyping import (
    make_joint_genotyping_jobs,
    JointGenotyperTool,
)
from cpg_pipes.jobs.seqr_loader import annotate_dataset_jobs, annotate_cohort_jobs
from cpg_pipes.jobs.somalier import check_pedigree_job
from cpg_pipes.jobs.vqsr import make_vqsr_jobs
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.types import CramPath, SequencingType

try:
    import utils
except ModuleNotFoundError:
    from . import utils


logger = logging.getLogger(__file__)


def _read_file(path: Path) -> str:
    with path.open() as f:
        return f.read().strip()


class TestJobs(unittest.TestCase):
    """
    Test individual jobs
    """

    @property
    def out_bucket(self):
        return utils.BASE_BUCKET / self.name / self.timestamp

    @property
    def tmp_bucket(self):
        return utils.BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def setUp(self):
        utils.setup_env()

        self.name = self._testMethodName
        self.timestamp = utils.timestamp()
        logger.info(f'Timestamp: {self.timestamp}')
        self.local_tmp_dir = tempfile.mkdtemp()

        self.pipeline = create_pipeline(
            name=self.name,
            description=self.name,
            analysis_dataset=utils.DATASET,
            namespace=Namespace.TEST,
        )
        self.sequencing_type = SequencingType.GENOME
        self.dataset = self.pipeline.create_dataset(utils.DATASET)
        sample_name = f'Test-{self.timestamp}'
        self.sample = self.dataset.add_sample(sample_name, sample_name)
        self.refs = self.pipeline.refs

        # Interval to take on chr20:
        self.chrom = 'chr20'
        self.locus1 = '5111495'
        self.locus2 = '5111607'
        self.interval = f'{self.chrom}-{self.locus1}-{self.locus2}'

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
        cram = self.pipeline.b.read_input_group(
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

    def test_align_fastq(self):
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
        self.pipeline.run(wait=True)

        self.assertTrue(
            (
                qc_bucket
                / 'duplicate-metrics'
                / f'{self.sample.id}-duplicate-metrics.csv'
            ).exists()
        )

        self.assertEqual(self.sample.id, _read_file(cram_details_paths['sample_name']))
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

    def test_align_cram(self):
        """
        Test alignment job on a CRAM input (tests realignment with bazam)
        """
        sid = 'CPG196519'
        cram_path = CramPath(utils.BASE_BUCKET / f'inputs/toy/cram/{sid}.cram')
        output_path = self.out_bucket / 'result.cram'
        jobs = align(
            self.pipeline.b,
            alignment_input=cram_path,
            output_path=output_path,
            sample_name=sid,
            refs=self.refs,
            realignment_shards_num=4,
        )
        cram_details_paths = self._job_get_cram_details(
            output_path,
            out_bucket=self.out_bucket / 'align',
            jobs=jobs,
        )
        self.pipeline.run(wait=True)

        self.assertEqual(sid, _read_file(cram_details_paths['sample_name']))
        cram_details = {k: _read_file(v) for k, v in cram_details_paths.items()}
        print(cram_details_paths)
        self.assertAlmostEqual(223007, int(cram_details['reads_num']), delta=50)
        self.assertAlmostEqual(
            222678,
            int(cram_details['reads_num_mapped_in_proper_pair']),
            delta=50,
        )

    def test_genotype_sample(self):
        """
        Test individual sample haplotype calling.
        """
        sid = 'CPG196519'
        cram_path = CramPath(
            utils.BASE_BUCKET / f'inputs/toy/cram_realigned/{sid}.cram'
        )
        out_gvcf_path = self.out_bucket / 'test.g.vcf.gz'
        jobs = produce_gvcf(
            self.pipeline.b,
            sample_name=sid,
            refs=self.refs,
            cram_path=cram_path,
            scatter_count=4,
            tmp_bucket=self.tmp_bucket,
            dragen_mode=True,
            output_path=out_gvcf_path,
            sequencing_type=self.sequencing_type,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        test_result_path = self._job_get_gvcf_header(out_gvcf_path, jobs)
        self.pipeline.run(wait=True)
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

        proj = self.pipeline.cohort.create_dataset(utils.DATASET)
        for sid in utils.SAMPLES:
            proj.add_sample(sid, sid)

        jobs = make_joint_genotyping_jobs(
            b=self.pipeline.b,
            out_vcf_path=out_vcf_path,
            out_siteonly_vcf_path=out_siteonly_vcf_path,
            refs=self.refs,
            gvcf_by_sid=utils.EXOME_1PCT_GVCF_BY_SID,
            tmp_bucket=self.tmp_bucket,
            overwrite=True,
            scatter_count=4,
            tool=JointGenotyperTool.GenotypeGVCFs,
            sequencing_type=self.pipeline.cohort.get_sequencing_type(),
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        test_result_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.pipeline.run(wait=True)
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
            b=self.pipeline.b,
            input_vcf_or_mt_path=siteonly_vcf_path,
            refs=self.refs,
            tmp_bucket=tmp_vqsr_bucket,
            gvcf_count=len(utils.SAMPLES),
            scatter_count=4,
            output_vcf_path=out_vcf_path,
            use_as_annotations=True,
            overwrite=True,
            sequencing_type=self.sequencing_type,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome5pct/calling_regions.interval_list',
        )
        res_path = self._job_get_gvcf_header(out_vcf_path, jobs)
        self.pipeline.run(wait=True)
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
            self.pipeline.b,
            vcf_path=siteonly_vcf_path,
            refs=self.refs,
            sequencing_type=self.sequencing_type,
            out_path=out_vcf_path,
            scatter_count=4,
            overwrite=False,
            hail_billing_project=utils.DATASET,
            hail_bucket=self.tmp_bucket,
            tmp_bucket=self.tmp_bucket,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )

        # Add test job
        test_j = self.pipeline.b.new_job('Parse VEP VCF')
        test_j.image(images.BCFTOOLS_IMAGE)
        test_j.command(
            f"""
        bcftools +split-vep {self.pipeline.b.read_input(str(out_vcf_path))} \
        -f '%CHROM:%POS %SYMBOL %BIOTYPE %MANE_SELECT %LoF %LoF_filter\n' \
        -i'BIOTYPE="protein_coding" & LoF_filter="ANC_ALLELE"' -s worst \
        > {test_j.output}
        """
        )
        test_j.depends_on(*jobs)
        test_out_path = self.out_bucket / f'{self.sample.id}.out'
        self.pipeline.b.write_output(test_j.output, str(test_out_path))

        # Run Batch
        self.pipeline.run(wait=True)

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
            self.pipeline.b,
            vcf=siteonly_vqsr_path,
            refs=self.refs,
            out_path=out_path,
            out_format='json',
        )

        # Add test job
        mane_transcript = 'NM_001009923.2'  # expected transcript on locus1
        test_j = self.pipeline.b.new_job('Parse VEP results')
        test_j.image(images.DRIVER_IMAGE)
        test_j.command(
            f"""
        cat {self.pipeline.b.read_input(str(out_path))} | zgrep {locus1} | \
        jq -r '.transcript_consequences[] | select(.mane_select=="{mane_transcript}") | .mane_select' \
        > {test_j.output}
        """
        )
        test_j.depends_on(j)
        test_out_path = self.out_bucket / f'{self.sample.id}.out'
        self.pipeline.b.write_output(test_j.output, str(test_out_path))

        # Run Batch
        self.pipeline.run(wait=True)

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
            b=self.pipeline.b,
            vep_results_paths=[vep_json_list_path],
            hail_billing_project=utils.DATASET,
            hail_bucket=self.tmp_bucket,
            out_path=out_path,
        )
        self.pipeline.run(wait=True)

        import hail as hl

        init_batch(utils.DATASET, self.tmp_bucket)
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
            self.pipeline.b,
            vcf_path=vcf_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_path,
            vep_ht_path=vep_ht_path,
            output_mt_path=out_mt_path,
            checkpoints_bucket=self.tmp_bucket / 'checkpoints',
            sequencing_type=self.sequencing_type,
            hail_billing_project=utils.DATASET,
            hail_bucket=self.tmp_bucket,
        )
        self.pipeline.run(wait=True)

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
            self.pipeline.b,
            vcf_path=vcf_path,
            siteonly_vqsr_vcf_path=siteonly_vqsr_path,
            vep_ht_path=vep_ht_path,
            output_mt_path=out_mt_path,
            checkpoints_bucket=self.tmp_bucket / 'checkpoints',
            sequencing_type=self.sequencing_type,
            hail_billing_project=utils.DATASET,
            hail_bucket=self.tmp_bucket,
        )
        self.pipeline.run(wait=True)

        # Testing
        import hail as hl

        init_batch(utils.DATASET, self.tmp_bucket)
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
            b=self.pipeline.b,
            mt_path=mt_path,
            sample_ids=utils.SAMPLES[:3],
            tmp_bucket=self.tmp_bucket,
            output_mt_path=out_mt_path,
            hail_billing_project=utils.DATASET,
            hail_bucket=self.tmp_bucket,
        )
        self.pipeline.run(wait=True)

        # Testing
        self.assertTrue(out_mt_path.exists())
        import hail as hl

        init_batch(utils.DATASET, self.tmp_bucket)
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
            self.pipeline.b,
            samples_file=self.pipeline.b.read_input(
                str(inputs_bucket / 'somalier-samples.tsv')
            ),
            pairs_file=self.pipeline.b.read_input(
                str(inputs_bucket / 'somalier-pairs.tsv')
            ),
        )
        self.pipeline.run(wait=True)
