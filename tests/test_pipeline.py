"""
Test building and running pipelines.
"""

import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch, Mock

from cpg_pipes import Namespace, Path, to_path
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.providers.cpg import CpgStorageProvider
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.types import SequencingType, CramPath

# Importing `seqr_loader` will make pipeline use all its stages by default.
from pipelines import seqr_loader
from pipelines.seqr_loader import AnnotateDataset

try:
    from .utils import BASE_BUCKET, SAMPLES, setup_env, DATASET, SUBSET_CRAM_BY_SID
except ImportError:
    from utils import BASE_BUCKET, SAMPLES, setup_env, DATASET, SUBSET_CRAM_BY_SID  # type: ignore


class UnittestStorageProvider(CpgStorageProvider):
    """
    Align and GenotypeSample stages write cram and gvcf into the datasets
    main bucket, without versioning. We need to override this behaviour to
    support writing into the test output directory.
    """

    def __init__(self, test_output_bucket: Path):
        super().__init__()
        self.test_output_bucket = test_output_bucket

    def _dataset_bucket(
        self,
        dataset: str,
        namespace: Namespace,
        suffix: str = None,
    ) -> Path:
        """
        Overiding main bucket name
        """
        path = self.test_output_bucket
        if suffix:
            path /= suffix
        return path


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

        self.realignment_shards_num = 4
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
        Mocking joint-calling outputs. Toy CRAM/GVCF don't produce enough variant data
        for AS-VQSR to work properly: gatk would throw a "Bad input: Values for
        AS_ReadPosRankSum annotation not detected for ANY training variant in the
        input callset", like in this Batch:
        https://batch.hail.populationgenomics.org.au/batches/46981/jobs/322.
        So instead, we are using mocks to plug in a larger file here.
        """
        jc_vcf = BASE_BUCKET / 'inputs/chr20/genotypegvcfs/joint-called.vcf.gz'
        siteonly_vcf = to_path(str(jc_vcf).replace('.vcf.gz', '-siteonly.vcf.gz'))
        from cpg_pipes.pipeline.pipeline import StageOutput

        original_fn = StageOutput.as_path

        def _mock_as_path(self_, id_=None):
            if self_.stage.name == 'JointGenotyping':
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
        first_stage: str | None = None,
        last_stage: str | None = None,
        realignment_shards_num: int = None,
        hc_intervals_num: int = None,
        jc_intervals_num: int = None,
        vep_intervals_num: int = None,
    ):
        setup_env()

        # Mocking elastic search password for the full dry run test:
        seqr_loader._read_es_password = Mock(return_value='TEST')

        pipeline = Pipeline(
            namespace=Namespace.TEST,
            name=self._testMethodName,
            analysis_dataset_name=DATASET,
            check_intermediates=False,
            storage_provider=UnittestStorageProvider(self.out_bucket),
            first_stage=first_stage,
            last_stage=last_stage,
            version=self.timestamp,
            config=dict(
                realignment_shards_num=realignment_shards_num
                or self.realignment_shards_num,
                hc_intervals_num=hc_intervals_num or self.hc_intervals_num,
                jc_intervals_num=jc_intervals_num or self.jc_intervals_num,
                vep_intervals_num=vep_intervals_num or self.vep_intervals_num,
            ),
        )
        self.datasets = [pipeline.create_dataset(DATASET)]
        for ds in self.datasets:
            for s_id in self.sample_ids:
                s = ds.add_sample(s_id, s_id)
                s.alignment_input = CramPath(SUBSET_CRAM_BY_SID[s.id])
                s.sequencing_type = seq_type
        return pipeline

    def test_dry(self):
        """
        Constucting a pipeline and submitting it to Hail Batch with dry_run.

        With dry_run, Hail Batch prints all code with a single print() call.
        Thus, we capture `builtins.print`, and verify that it has the expected
        job commands passed to it.

        """
        pipeline = self._setup_pipeline()

        with patch('builtins.print') as mock_print:
            result = pipeline.run(dry_run=True)
            self.assertEqual('success', result.status()['state'])

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
        self.assertEqual(
            _cnt('HaplotypeCaller'), len(self.sample_ids) * self.hc_intervals_num
        )
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

    def test_joint_calling(self):
        """
        Stages up to joint calling.
        """
        pipeline = self._setup_pipeline(
            last_stage=JointGenotyping.__name__,
            seq_type=SequencingType.TOY,
        )
        result = pipeline.run(dry_run=False, wait=True)
        self.assertEqual('success', result.status()['state'])

    def test_after_joint_calling(self):
        """
        Stages after joint-calling (running separately, because
        VQSR needs more inputs than provided by toy regions)
        """
        pipeline = self._setup_pipeline(
            first_stage=Vqsr.__name__,
            last_stage=AnnotateDataset.__name__,
            seq_type=SequencingType.WGS,
        )
        # Mocking joint-calling outputs. Toy CRAM/GVCF don't produce enough variant
        # data for AS-VQSR to work properly: gatk would throw a "Bad input: Values for
        # AS_ReadPosRankSum annotation not detected for ANY training variant in the
        # input callset", like in this Batch:
        # https://batch.hail.populationgenomics.org.au/batches/46981/jobs/322.
        # So instead, we are copying a larger file into JointGenotyping expected output.
        jc_vcf = BASE_BUCKET / 'inputs/chr20/genotypegvcfs/joint-called.vcf.gz'
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
