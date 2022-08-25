"""
Test Hail Query functions.
"""

import logging
import shutil
import tempfile
import unittest
from pprint import pprint

import hail as hl
from cpg_utils import to_path
from cpg_utils.workflows.utils import timestamp
from cpg_utils.hail_batch import init_batch

from jobs.query_scripts.seqr_loader import (
    annotate_cohort,
    subset_mt_to_samples,
    annotate_dataset_mt,
    load_vqsr,
)
from jobs.query_scripts.vep import vep_json_to_ht

logger = logging.getLogger(__file__)


try:
    import utils
except ModuleNotFoundError:
    from . import utils


class TestQuery(unittest.TestCase):
    """
    Test Hail Query functions.
    """

    @property
    def out_bucket(self):
        """Property to allow re-setting the timestamp in a method"""
        return utils.BASE_BUCKET / self.name / self.timestamp

    @property
    def tmp_bucket(self):
        """Property to allow re-setting the timestamp in a method"""
        return utils.BASE_BUCKET / 'tmp' / self.name / self.timestamp

    def setUp(self):
        """Called for each test method"""
        self.name = self._testMethodName
        self.timestamp = timestamp()
        self.local_tmp_dir = tempfile.mkdtemp()
        self.sequencing_type = 'genome'
        init_batch()
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

    def test_export_vqsr_vcf(self):
        """
        Test exporting AS-VQSR site-only VCF (it's done test_seqr_loader_annotate_cohort),
        but we want to test it separately.
        """
        siteonly_vqsr_path = to_path(
            'gs://cpg-fewgenomes-test/unittest/inputs/vqsr_parse/siteonly-regions.vcf.gz'
        )
        ht = load_vqsr(str(siteonly_vqsr_path))
        results = ht.info.collect()
        pprint(results)
        self.assertEqual(results[0].AS_FilterStatus, 'PASS')
        self.assertEqual(results[1].AS_FilterStatus, 'VQSRTrancheINDEL99.50to99.90')
        # Checking that AS-VQSR-used INFO/AS_* fields are parsed correctly
        self.assertEqual(results[0].AS_ReadPosRankSum, -0.3)
        self.assertEqual(results[0].AS_MQRankSum, -3.8)
        self.assertEqual(results[1].AS_MQ, 11.37)
        self.assertEqual(results[1].AS_FS, 3.592)
        self.assertEqual(results[1].AS_SOR, 1.267)

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
            f'gs://cpg-fewgenomes-test/unittest/inputs/chr20/vep/{self.interval}.ht'
        )
        out_mt_path = self.out_bucket / 'seqr_loader' / f'cohort-{self.interval}.mt'
        annotate_cohort(
            str(vcf_path),
            str(siteonly_vqsr_path),
            vep_ht_path=str(vep_ht_path),
            out_mt_path=str(out_mt_path),
            checkpoints_bucket=str(self.tmp_bucket / 'seqr_loader' / 'checkpoints'),
            sequencing_type=self.sequencing_type.upper(),
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
        subset_mt_path = self.tmp_bucket / 'seqr_loader' / f'subset-{self.interval}.mt'
        out_mt_path = self.out_bucket / 'seqr_loader' / f'dataset-{self.interval}.mt'
        subset_mt_to_samples(
            str(mt_path),
            utils.SAMPLES[:3],
            str(subset_mt_path),
        )
        annotate_dataset_mt(
            str(subset_mt_path), str(out_mt_path), str(self.tmp_bucket / 'seqr_loader')
        )
        # Testing
        mt = hl.read_matrix_table(str(out_mt_path))
        mt.rows().show()
        self.assertListEqual(mt.s.collect(), utils.SAMPLES[:3])
        self.assertSetEqual(
            set(mt.samples_gq['20_to_25'].collect()[0]), {'CPG196519', 'CPG196527'}
        )
        self.assertSetEqual(set(mt.samples_ab['40_to_45'].collect()[0]), {'CPG196535'})
