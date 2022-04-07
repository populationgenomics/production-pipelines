"""
Test building and running pipelines.
"""

import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch, Mock

from cpg_pipes import Namespace, to_path
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.refdata import RefData
from cpg_pipes.types import SequencingType

try:
    from .utils import setup_env, BASE_BUCKET, DATASET, SAMPLES
except ImportError:
    from utils import setup_env, BASE_BUCKET, DATASET, SAMPLES  # type: ignore


class TestPipeline(unittest.TestCase):
    """
    Test building and running pipelines.
    """
    def setUp(self):
        """
        Setting parameters, creating local tmp dir.
        """
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.out_bucket = BASE_BUCKET / self.name / self.timestamp
        self.tmp_bucket = self.out_bucket / 'tmp'
        self.local_tmp_dir = to_path(tempfile.mkdtemp())
        self.sample_ids = SAMPLES[:3]
        
        self.hc_intervals_num = 4
        self.jc_intervals_num = 4
        self.vep_intervals_num = 4
        
    def tearDown(self) -> None:
        """
        Removing local tmp dir
        """
        shutil.rmtree(self.local_tmp_dir)
        
    @staticmethod
    def _mock_joint_calling():
        """
        Mocking joint-calling outputs. Toy FASTQs don't produce enough variant data
        for AS-VQSR to work properly: gatk would throw a "Bad input: Values for 
        AS_ReadPosRankSum annotation not detected for ANY training variant in the 
        input callset", like in this Batch:
        https://batch.hail.populationgenomics.org.au/batches/46981/jobs/322.
        So instead, we are using mocks to plug in a larger file here.
        """
        jc_vcf = to_path(
            'gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/'
            'joint-called.vcf.gz'
        )
        siteonly_vcf = to_path(str(jc_vcf).replace('.vcf.gz', '-siteonly.vcf.gz'))
        from cpg_pipes.pipeline.pipeline import StageOutput
        original_fn = StageOutput.as_path

        def _mock_as_path(self_, id_=None):
            if self_.stage.name == 'JointGenotypingStage':
                d = {
                    'vcf': jc_vcf, 
                    'siteonly': siteonly_vcf,
                }
                return d[id_] if id_ else d
            return original_fn(self_, id_)

        StageOutput.as_path = _mock_as_path
        
    def _setup_pipeline(
        self, 
        seq_type=SequencingType.WGS,
        last_stage=None,
    ):
        setup_env()
        from cpg_pipes import benchmark
        # Use the seqr_loader stages. Importing it will make sure all its stages
        # are used by default:
        from pipelines import seqr_loader  # noqa: F401

        # Mocking elastic search password for the full dry run test:
        seqr_loader._read_es_password = Mock(return_value='TEST')
        
        # Mocking joint-calling results:
        self._mock_joint_calling()

        pipeline = create_pipeline(
            name=self._testMethodName,
            description=self._testMethodName,
            analysis_dataset=DATASET,
            namespace=Namespace.TEST,
            check_intermediates=False,
            check_expected_outputs=False,
            last_stage=last_stage,
            config=dict(
                hc_intervals_num=self.hc_intervals_num,
                jc_intervals_num=self.jc_intervals_num,
                vep_intervals_num=self.vep_intervals_num,
            )
        )
        self.datasets = [pipeline.add_dataset(DATASET)]
        for ds in self.datasets:
            for s_id in self.sample_ids:
                s = ds.add_sample(s_id, s_id)
                s.alignment_input = benchmark.tiny_fq
                s.sequencing_type = seq_type
        return pipeline

    def test_dry(self):
        """
        Constucting a pipeline and submitting it to Hail Batch with dry_run.

        With dry_run, Hail Batch prints all code with a single print() call.
        Thus, we capture `builtins.print`, and verify that it has the expected
        job commands passed to it.

        # Mocking all hail methods (hail.hadoop_open in this case) so we don't 
        # have to initialize hail, which steals a few seconds from this test
        # which is supposed to be quick.
        # """
        # with patch('cpg_pipes.utils.exists') as mock_exists:
        #     mock_exists.return_value = True
        pipeline = self._setup_pipeline()

        with patch('builtins.print') as mock_print:
            pipeline.run(dry_run=True)
        
            # print() should be called only once:
            self.assertEqual(1, mock_print.call_count)
    
            # all job commands would be contained in this one print call:
            out = mock_print.call_args_list[0][0][0]

        sys.stdout.write(out)
        lines = out.split('\n')
        
        def _cnt(item: str) -> int:
            """Number of lines that start with item"""
            return len([line for line in lines if line.strip().startswith(item)])

        self.assertEqual(_cnt('bwa mem'), len(self.sample_ids) * 2)
        self.assertEqual(_cnt('HaplotypeCaller'), len(self.sample_ids) * self.hc_intervals_num)
        self.assertEqual(_cnt('ReblockGVCF'), len(self.sample_ids))
        self.assertEqual(_cnt('GenotypeGVCFs'), self.jc_intervals_num)
        self.assertEqual(_cnt('GenomicsDBImport'), self.jc_intervals_num)
        self.assertEqual(_cnt('MakeSitesOnlyVcf'), self.jc_intervals_num)
        self.assertEqual(_cnt('VariantRecalibrator'), 2)
        self.assertEqual(_cnt('ApplyVQSR'), 2)
        self.assertEqual(_cnt('GatherVcfsCloud'), 2)
        self.assertEqual(_cnt('vep '), self.vep_intervals_num)
        self.assertEqual(_cnt('vep_json_to_ht('), 1)
        self.assertEqual(_cnt('annotate_cohort('), 1)
        self.assertEqual(_cnt('annotate_dataset_mt('), len(self.datasets))
        self.assertEqual(_cnt('hailctl dataproc submit'), 1)

    def test_wgs(self):
        """
        WGS seqr-loading pipeline.
        """
        pipeline = self._setup_pipeline(
            last_stage='AnnotateDataset'
        )
        pipeline.run(dry_run=False, wait=True)

    def test_exome(self):
        """
        Exome seqr-loading pipeline.
        """
        pipeline = self._setup_pipeline(
            last_stage='AnnotateDataset',
            seq_type=SequencingType.EXOME,
        )
        pipeline.run(dry_run=False, wait=True)


if __name__ == '__main__':
    unittest.main()
