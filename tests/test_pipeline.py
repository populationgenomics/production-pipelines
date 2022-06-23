"""
Test building and running pipelines.
"""
import shutil
import sys
import tempfile
import unittest
from unittest import skip
from unittest.mock import patch, Mock

from cpg_utils.config import get_config, set_config_paths

from cpg_pipes import to_path, get_package_path
from cpg_pipes.pipeline.pipeline import Pipeline, Stage, StageDecorator
from cpg_pipes.stages.genotype_sample import GenotypeSample
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.types import CramPath
from cpg_pipes.stages import seqr_loader
from cpg_pipes.utils import timestamp

try:
    import utils
except ModuleNotFoundError:
    from . import utils


class TestPipeline(unittest.TestCase):
    """
    Test building and running pipelines.
    """

    def setUp(self):
        """
        Setting parameters, creating local tmp dir.
        """
        self.name = self._testMethodName
        self.timestamp = timestamp()
        self.out_bucket = utils.BASE_BUCKET / self.name / self.timestamp
        self.tmp_bucket = self.out_bucket / 'tmp'
        self.local_tmp_dir = tempfile.mkdtemp()
        self.sample_ids = utils.SAMPLES[:3]
        
        self.config_path = self.tmp_bucket / 'config.toml'

    def setup_env(self, intervals_path: str | None = None):
        """
        Mock analysis-runner environment.
        """
        conf = {'workflow': {'version': self.timestamp}}
        if intervals_path:
            conf['workflow']['intervals_path'] = intervals_path
        utils.setup_env(self.tmp_bucket, conf)

    def tearDown(self) -> None:
        """
        Removing local tmp dir
        """
        shutil.rmtree(self.local_tmp_dir)

    def _setup_pipeline(self, stages: list[StageDecorator]):
        # Mocking elastic search password for the full dry run test:
        seqr_loader.es_password = Mock(return_value='TEST')
        pipeline = Pipeline(name=self._testMethodName, stages=stages)
        self.datasets = [pipeline.create_dataset(utils.DATASET)]
        for ds in self.datasets:
            for s_id in self.sample_ids:
                s = ds.add_sample(s_id, s_id)
                s.alignment_input_by_seq_type[utils.SEQ_TYPE] = \
                    CramPath(utils.TOY_CRAM_BY_SID[s.id])
        return pipeline

    def test_dry(self):
        """
        Constucting a pipeline and submitting it to Hail Batch with dry_run.

        With dry_run, Hail Batch prints all code with a single print() call.
        Thus, we capture `builtins.print`, and verify that it has the expected
        job commands passed to it.
        """
        self.setup_env()
        pipeline = self._setup_pipeline(stages=[seqr_loader.LoadToEs])

        with patch('builtins.print') as mock_print:
            with patch.object(Stage, '_outputs_are_reusable') as mock_reusable:
                mock_reusable.return_value = False
                pipeline.run(dry_run=True, force_all_implicit_stages=True)

            # print() should be called only once:
            self.assertEqual(mock_print.call_count, 1)

            # all job commands would be contained in this one print call:
            out = mock_print.call_args_list[0][0][0]

        sys.stdout.write(out)
        lines = out.split('\n')

        def _cnt(item: str) -> int:
            """Number of lines that start with item"""
            return len([line for line in lines if line.strip().startswith(item)])

        self.assertEqual(
            _cnt('dragen-os'), len(self.sample_ids) * get_config()['workflow']['realignment_shards_num']
        )
        self.assertEqual(
            _cnt('HaplotypeCaller'), len(self.sample_ids) * get_config()['workflow']['hc_intervals_num']
        )
        self.assertEqual(_cnt('ReblockGVCF'), len(self.sample_ids))
        self.assertEqual(_cnt('GenotypeGVCFs'), get_config()['workflow']['jc_intervals_num'])
        self.assertEqual(_cnt('GenomicsDBImport'), get_config()['workflow']['jc_intervals_num'])
        self.assertEqual(_cnt('MakeSitesOnlyVcf'), get_config()['workflow']['jc_intervals_num'])
        # Indel + SNP create model + SNP scattered
        self.assertEqual(_cnt('VariantRecalibrator'), 2 + get_config()['workflow']['jc_intervals_num'])
        # Twice to each interva: apply indels, apply SNPs
        self.assertEqual(_cnt('ApplyVQSR'), 2 * get_config()['workflow']['jc_intervals_num'])
        # Gather JC, gather siteonly JC, gather VQSR
        self.assertEqual(_cnt('GatherVcfsCloud'), 3)
        self.assertEqual(_cnt('vep '), get_config()['workflow']['vep_intervals_num'])
        self.assertEqual(_cnt('vep_json_to_ht('), 1)
        self.assertEqual(_cnt('annotate_cohort('), 1)
        self.assertEqual(_cnt('annotate_dataset_mt('), len(self.datasets))
        self.assertEqual(_cnt('hailctl dataproc submit'), 1)

    @skip('Running only dry tests in ths module')
    def test_joint_calling(self):
        """
        Stages up to joint calling.
        """
        self.setup_env(
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        pipeline = self._setup_pipeline(
            stages=[GenotypeSample.__name__, JointGenotyping.__name__],
        )
        result = pipeline.run(dry_run=False, wait=True)
        self.assertEqual('success', result.status()['state'])

    @skip('Running only dry tests in ths module. VEP takes too long')
    def test_after_joint_calling(self):
        """
        VQSR needs more inputs than provided by toy regions.
        """
        self.setup_env(
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome5pct/calling_regions.interval_list',
        )
        pipeline = self._setup_pipeline(
            stages=[Vqsr.__name__, seqr_loader.AnnotateDataset.__name__],
        )
        # Mocking joint-calling outputs. Toy CRAM/GVCF don't produce enough variant
        # data for AS-VQSR to work properly: gatk would throw a "Bad input: Values for
        # AS_ReadPosRankSum annotation not detected for ANY training variant in the
        # input callset", like in this Batch:
        # https://batch.hail.populationgenomics.org.au/batches/46981/jobs/322.
        # So instead, we are copying a larger file into JointGenotyping expected output.
        jc_vcf = utils.BASE_BUCKET / 'inputs/exome5pct/9samples-joint-called.vcf.gz'

        siteonly_vcf = to_path(str(jc_vcf).replace('.vcf.gz', '-siteonly.vcf.gz'))
        expected_output = JointGenotyping(pipeline).expected_outputs(pipeline.cohort)
        jc_vcf.copy(expected_output['vcf'], force_overwrite_to_cloud=True)
        to_path(str(jc_vcf) + '.tbi').copy(
            str(expected_output['vcf']) + '.tbi', force_overwrite_to_cloud=True
        )
        siteonly_vcf.copy(expected_output['siteonly'], force_overwrite_to_cloud=True)
        to_path(str(siteonly_vcf) + '.tbi').copy(
            str(expected_output['siteonly']) + '.tbi', force_overwrite_to_cloud=True
        )

        result = pipeline.run(dry_run=False, wait=True)

        self.assertEqual('success', result.status()['state'])


if __name__ == '__main__':
    unittest.main()
