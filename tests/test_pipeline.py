import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

try:
    from .utils import setup_env, BASE_BUCKET, PROJECT, SAMPLES
except ImportError:
    from utils import setup_env, BASE_BUCKET, PROJECT, SAMPLES  # type: ignore


class TestPipeline(unittest.TestCase):
    """
    Test Pipeline class
    """
    def setUp(self):
        self.name = self._testMethodName
        self.timestamp = time.strftime('%Y%m%d-%H%M')
        self.out_bucket = f'{BASE_BUCKET}/{self.name}/{self.timestamp}'
        self.tmp_bucket = f'{self.out_bucket}/tmp'
        self.local_tmp_dir = tempfile.mkdtemp()
        self.sample_name = f'Test-{self.timestamp}'
        self.sample_ids = SAMPLES[:3]

    def tearDown(self) -> None:
        shutil.rmtree(self.local_tmp_dir)

    @patch('hail.hadoop_open', lambda a, b: tempfile.TemporaryFile('w'))
    def test_pipeline(self):
        """
        Constucting a pipeline and submitting it to Hail Batch with dry_run.

        With dry_run, Hail Batch prints all code with a single print() call.
        Thus we capture builtins.print, and verify that it has the expected
        job commands passed to it.
        
        Mocking all hail methods (hail.hadoop_open in this case) so we don't 
        have to initialize hail, which steals a few seconds from this test
        which is supposed to be quick.
        """
        setup_env()
        from cpg_pipes import benchmark
        from cpg_pipes import ref_data
        from cpg_pipes.pipeline.pipeline import Pipeline
        # Use the seqr_loader stages. Importing it will make sure all its stages
        # are used by default:
        from pipelines import seqr_loader  # noqa: F401

        pipeline = Pipeline(
            name=self._testMethodName,
            description=self._testMethodName,
            analysis_project=PROJECT,
            output_version='v0',
            namespace='test',
            check_smdb_seq=False,
            config=dict(output_projects=[PROJECT]),
            dry_run=True,
        )
        project = pipeline.cohort.add_project(PROJECT)
        for s_id in self.sample_ids:
            s = project.add_sample(s_id, s_id)
            s.alignment_input = benchmark.tiny_fq

        with patch('builtins.print') as mock_print:
            pipeline.submit_batch()
        
            # print() should be called only once:
            self.assertEqual(1, mock_print.call_count)
    
            # all job commands would be contained in this one print call:
            out = mock_print.call_args_list[0][0][0]

        sys.stdout.write(out)
        lines = out.split('\n')
        
        def _cnt(item: str) -> int:
            """Number of lines that start with item"""
            return len([line for line in lines if line.strip().startswith(item)])

        self.assertEqual(_cnt('bwa mem'), len(self.sample_ids))
        self.assertEqual(_cnt('HaplotypeCaller'), len(self.sample_ids))
        self.assertEqual(_cnt('ReblockGVCF'), len(self.sample_ids))
        self.assertEqual(_cnt('GenotypeGVCFs'), ref_data.NUMBER_OF_GENOMICS_DB_INTERVALS)
        self.assertEqual(_cnt('GenomicsDBImport'), ref_data.NUMBER_OF_GENOMICS_DB_INTERVALS)
        self.assertEqual(_cnt('MakeSitesOnlyVcf'), ref_data.NUMBER_OF_GENOMICS_DB_INTERVALS)
        self.assertEqual(_cnt('VariantRecalibrator'), 2)
        self.assertEqual(_cnt('ApplyVQSR'), 2)
        self.assertEqual(_cnt('GatherVcfsCloud'), 2)
        self.assertEqual(
            _cnt('hailctl dataproc submit'), 
            1 + len(pipeline.config['output_projects']) * 3
        )
