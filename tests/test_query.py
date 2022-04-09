"""
Test Hail Query functions.
"""

import shutil
import tempfile
import time
import unittest
from unittest import skip

import hail as hl

from cpg_pipes import to_path, Namespace
from cpg_pipes import hailquery
from cpg_pipes.providers.cpg import CpgStorageProvider
from cpg_pipes.query.seqr_loader import annotate_cohort, subset_mt_to_samples, \
    annotate_dataset_mt
from cpg_pipes.query.vep import vep_json_to_ht
from cpg_pipes.refdata import RefData
from cpg_pipes.types import SequencingType, logger

try:
    from .utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID
except ImportError:
    from utils import BASE_BUCKET, DATASET, SAMPLES, SUBSET_GVCF_BY_SID  # type: ignore


class TestQuery(unittest.TestCase):
    """
    Test Hail Query functions.
    """
    @property
    def out_bucket(self):
        """Property to allow re-setting the timestamp in a method"""
        return BASE_BUCKET / self.name / self.timestamp

    @property
    def tmp_bucket(self):
        """Property to allow re-setting the timestamp in a method"""
        return BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def setUp(self):
        """Called for each test method"""
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.local_tmp_dir = tempfile.mkdtemp()
        self.sequencing_type = SequencingType.WGS
        self.refs = RefData(CpgStorageProvider().get_ref_bucket())
        hailquery.init_batch(DATASET, self.tmp_bucket)
        # Interval to take on chr20:
        self.chrom = 'chr20'
        self.locus1 = '5111495'
        self.locus2 = '5111607'
        self.interval = f'{self.chrom}-{self.locus1}-{self.locus2}'
        
    def tearDown(self) -> None:
        """Remove tmp dir"""
        shutil.rmtree(self.local_tmp_dir)

    def test_vep_json_to_ht(self):
        """
        Test parsing VEP JSON into a hail table.
        """
        vep_json_list_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/vep/'
            f'{self.interval}.json_list'
        )
        out_path = self.out_bucket / 'vep_ht' / f'{self.interval}.ht'

        vep_json_to_ht(
            vep_results_paths=[str(vep_json_list_path)],
            out_path=str(out_path),
        )
        logger.info(f'Written VEP annotations into {out_path}')

        # Testing
        ht = hl.read_table(str(out_path))
        interval = hl.parse_locus_interval(
            f'{self.chrom}:{self.locus1}-{int(self.locus1) + 1}'
        )
        ht = hl.filter_intervals(ht, [interval])
        transcript = ht.vep.transcript_consequences.mane_select.collect()[0][1]
        lof = ht.vep.transcript_consequences.lof.collect()[0][1]
        self.assertEqual(lof, 'LC')
        self.assertEqual(transcript, 'NM_001009923.2')

    def test_seqr_loader_annotate_cohort(self):
        """
        Test seqr loader annotation: cohort VCF into a mt.
        """
        vcf_path = (
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'joint-called-{self.interval}.vcf.gz'
        )
        siteonly_vqsr_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'siteonly-vqsr-{self.interval}.vcf.gz'
        )
        vep_ht_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'vep/{self.interval}.ht'
        )
        out_mt_path = (
            self.out_bucket / 'seqr_loader' / f'cohort-{self.interval}.mt'
        )
        annotate_cohort(
            str(vcf_path),
            str(siteonly_vqsr_path),
            vep_ht_path=str(vep_ht_path),
            out_mt_path=str(out_mt_path),
            checkpoints_bucket=str(self.tmp_bucket / 'seqr_loader' / 'checkpoints'),
            sequencing_type=self.sequencing_type.value.upper(),
        )
        # Testing
        mt = hl.read_matrix_table(str(out_mt_path))
        mt.rows().show()
        self.assertListEqual(mt.topmed.AC.collect(), [20555, 359, 20187])
        self.assertSetEqual(set(mt.geneIds.collect()[0]), {'ENSG00000089063'})

    def test_seqr_loader_annotate_dataset(self):
        mt_path = to_path(
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/'
            f'seqr_loader/cohort-{self.interval}.mt'
        )
        subset_mt_path = (
            self.tmp_bucket / 'seqr_loader' / f'subset-{self.interval}.mt'
        )
        out_mt_path = (
            self.out_bucket / 'seqr_loader' / f'dataset-{self.interval}.mt'
        )
        subset_mt_to_samples(
            str(mt_path),
            SAMPLES[:3],
            str(subset_mt_path),
        )
        annotate_dataset_mt(
            str(subset_mt_path),
            str(out_mt_path),
            str(self.tmp_bucket / 'seqr_loader')
        )
        # Testing
        mt = hl.read_matrix_table(str(out_mt_path))
        mt.rows().show()
        self.assertListEqual(mt.s.collect(), SAMPLES[:3])
        self.assertSetEqual(
            set(mt.samples_gq['20_to_25'].collect()[0]), {'CPG196519', 'CPG196527'}
        )
        self.assertSetEqual(
            set(mt.samples_ab['40_to_45'].collect()[0]), {'CPG196535'}
        )

    @skip('Not implemented in Batch backend')
    def test_vcf_combiner(self):
        from cpg_pipes.targets import Cohort
        dataset = Cohort(
            analysis_dataset_name=DATASET,
            namespace=Namespace.TEST,
            storage_provider=CpgStorageProvider()
        ).add_dataset(DATASET)
        for sid in SAMPLES:
            dataset.add_sample(sid)
    
        out_mt_path = self.out_bucket / 'combined.mt'
        
        hl.experimental.run_combiner(
            [str(s.get_gvcf_path().path) for s in dataset.get_samples()],
            sample_names=[s.id for s in dataset.get_samples()],
            out_file=str(out_mt_path),
            reference_genome=RefData.genome_build,
            use_genome_default_intervals=True,
            tmp_path=self.tmp_bucket,
            overwrite=True,
            key_by_locus_and_alleles=True,
        )

        mt = hl.read_matrix_table(out_mt_path)
        logger.info(
            f'Written {mt.cols().count()} samples to {out_mt_path}, '
            f'n_partitions={mt.n_partitions()}'
        )
        self.assertSetEqual(set(mt.s.collect()), SAMPLES)
