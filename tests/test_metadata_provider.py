"""
Testing TSV metadata provider. TSV file must contain one required column "sample".
"""
import unittest
from io import StringIO

from cpg_pipes.pipeline.dataset import Cohort, Sex
from cpg_pipes.pipeline.metadata_provider import CsvMetadataProvider
from cpg_pipes.storage import Namespace
from cpg_pipes.cpg.storage import CPGStorageProvider


DATASET = 'fewgenomes'
SAMPLE1 = 'CPG01'
SAMPLE2 = 'CPG02'
SAMPLE3 = 'CPG03'
EXTID1 = 'MYSAMPLE1'
EXTID2 = 'MYSAMPLE2'
EXTID3 = 'MYSAMPLE3'
FQ = 'gs://cpg-fewgenomes-test-upload/MYSAMPLE1_M00{}_R{}.fastq.gz'
CRAM = 'gs://cpg-fewgenomes-test-upload/MYSAMPLE2.cram'


class TestTsvProvider(unittest.TestCase):
    """
    Test TSV metadat provider.
    """

    def test_csv_provider(self):
        tsv_contents = f"""\
dataset,sample,external_id,fqs_r1,fqs_r2,cram,sex
{DATASET},{SAMPLE1},{EXTID1},{FQ.format(1, 1)}|{FQ.format(2, 1)},{FQ.format(2, 1)}|{FQ.format(2, 2)},,M
{DATASET},{SAMPLE2},{EXTID2},,,{CRAM},
{DATASET},{SAMPLE3},{EXTID3},,,,
        """.strip()

        with StringIO(tsv_contents) as fp:
            provider = CsvMetadataProvider(fp)
            
        cohort = Cohort(
            name='test',
            analysis_dataset_name=DATASET,
            namespace=Namespace.TEST,
            storage_provider=CPGStorageProvider(),
        )
        provider.populate_cohort(
            cohort, 
            do_check_seq_existence=False,
            skip_samples=[EXTID3],
        )
        
        dss = cohort.get_datasets()
        self.assertEqual(len(dss), 1)
        ds = dss[0]
        self.assertEqual(ds.name, f'{DATASET}-test')
        ss = cohort.get_datasets()[0].get_samples()
        self.assertEqual(len(ss), 2)
        s1 = ss[0]
        self.assertEqual(s1.external_id, EXTID1)
        self.assertEqual(len(s1.alignment_input), 2)
        self.assertEqual(s1.pedigree.sex, Sex.MALE)
        s2 = ss[1]
        self.assertEqual(s2.external_id, EXTID2)
        self.assertTrue(s2.alignment_input.ext == 'cram')
        self.assertEqual(s2.pedigree.sex, Sex.UNKNOWN)


if __name__ == '__main__':
    unittest.main()
