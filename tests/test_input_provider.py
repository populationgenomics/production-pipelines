"""
Testing input providers.
"""

import unittest
from io import StringIO
from unittest import skip

from cpg_pipes import Namespace
from cpg_pipes.providers.cpg.smdb import SMDB
from cpg_pipes.targets import Cohort, Sex
from cpg_pipes.providers.inputs import CsvInputProvider
from cpg_pipes.providers.cpg import CpgStorageProvider, SmdbInputProvider
from cpg_pipes.types import SequencingType


class TestInputProvider(unittest.TestCase):
    """
    Testing inputs providers.
    """

    def test_csv_provider(self):
        """
        Test CSV/CSV input metadata provider. A file must contain one required
        column "sample".
        """
        dataset = 'fewgenomes'
        sample1 = 'CPG01'
        sample2 = 'CPG02'
        sample3 = 'CPG03'
        extid1 = 'MYSAMPLE1'
        extid2 = 'MYSAMPLE2'
        extid3 = 'MYSAMPLE3'
        fq = 'gs://cpg-fewgenomes-test-upload/MYSAMPLE1_M00{}_R{}.fastq.gz'
        cram = 'gs://cpg-fewgenomes-test-upload/MYSAMPLE2.cram'

        tsv_contents = f"""\
dataset,sample,external_id,fqs_r1,fqs_r2,cram,sex,seq_type
{dataset},{sample1},{extid1},{fq.format(1, 1)}|{fq.format(2, 1)},{fq.format(2, 1)}|{fq.format(2, 2)},,M,wgs
{dataset},{sample2},{extid2},,,{cram},,exome
{dataset},{sample3},{extid3},,,,,wgs
        """.strip()

        with StringIO(tsv_contents) as fp:
            provider = CsvInputProvider(fp, check_files=False)

        cohort = Cohort(
            name='test',
            analysis_dataset_name=dataset,
            namespace=Namespace.TEST,
            storage_provider=CpgStorageProvider(),
        )
        provider.populate_cohort(
            cohort,
            skip_samples=[extid3],
        )
        self.assertEqual(len(cohort.get_datasets()), 1)
        ds = cohort.get_datasets()[0]
        self.assertEqual(ds.name, f'{dataset}-test')
        self.assertEqual(len(ds.get_samples()), 2)
        s1 = ds.get_samples()[0]
        self.assertEqual(s1.external_id, extid1)
        self.assertEqual(len(s1.alignment_input), 2)
        self.assertEqual(s1.pedigree.sex, Sex.MALE)
        self.assertEqual(s1.sequencing_type, SequencingType.WGS)
        s2 = ds.get_samples()[1]
        self.assertEqual(s2.external_id, extid2)
        self.assertTrue(s2.alignment_input.ext == 'cram')
        self.assertEqual(s2.pedigree.sex, Sex.UNKNOWN)
        self.assertEqual(s2.sequencing_type, SequencingType.EXOME)

    @skip('Figure out SMDB permissions from GitHub workflows')
    def test_smdb_provider(self):
        """
        Test sample-metadata input provider.
        """
        dataset = 'acute-care'

        input_provider = SmdbInputProvider(SMDB())
        cohort = input_provider.populate_cohort(
            cohort=Cohort(
                analysis_dataset_name=dataset,
                namespace=Namespace.TEST,
                storage_provider=CpgStorageProvider(),
            ),
            dataset_names=[dataset],
        )

        self.assertEqual(len(cohort.get_datasets()), 1)
        ds = cohort.get_datasets()[0]
        self.assertEqual(ds.name, f'{dataset}-test')
        self.assertEqual(len(ds.get_samples()), 9)
        s1 = ds.get_samples()[0]
        s2 = ds.get_samples()[1]
        self.assertEqual(s1.external_id, '20W002328')
        self.assertEqual(len(s1.alignment_input), 2)
        self.assertEqual(s1.pedigree.sex, Sex.FEMALE)
        self.assertEqual(s1.pedigree.fam_id, s2.pedigree.fam_id)
        self.assertEqual(s1.pedigree.mom.id, s2.id)


if __name__ == '__main__':
    unittest.main()
