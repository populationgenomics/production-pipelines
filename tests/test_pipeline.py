"""
Test building and running pipelines.
"""

import shutil
import sys
import tempfile
import unittest
from unittest import skip
from unittest.mock import patch, Mock

from cpg_pipes import Namespace, Path, to_path
from cpg_pipes.pipeline.pipeline import Pipeline, Stage
from cpg_pipes.providers.cpg.images import CpgImages
from cpg_pipes.providers.cpg.refdata import CpgRefData
from cpg_pipes.providers.cpg.storage import CpgStorageProvider
from cpg_pipes.stages.genotype_sample import GenotypeSample
from cpg_pipes.stages.joint_genotyping import JointGenotyping
from cpg_pipes.stages.vqsr import Vqsr
from cpg_pipes.types import CramPath, AlignmentInput

# Importing `seqr_loader` will make pipeline use all its stages by default.
from pipelines import seqr_loader
from pipelines.seqr_loader import AnnotateDataset

try:
    import utils
except ModuleNotFoundError:
    from . import utils


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
        self.timestamp = utils.timestamp()
        self.out_bucket = utils.BASE_BUCKET / self.name / self.timestamp
        self.tmp_bucket = self.out_bucket / 'tmp'
        self.local_tmp_dir = to_path(tempfile.mkdtemp())
        self.sample_ids = utils.SAMPLES[:3]

        self.realignment_shards_num = 4
        self.hc_intervals_num = 4
        self.jc_intervals_num = 4
        self.vep_intervals_num = 4

    def tearDown(self) -> None:
        """
        Removing local tmp dir
        """
        shutil.rmtree(self.local_tmp_dir)

    def _setup_pipeline(
        self,
        first_stage: str | None = None,
        last_stage: str | None = None,
        realignment_shards_num: int = None,
        intervals_path: Path | None = None,
        hc_intervals_num: int = None,
        jc_intervals_num: int = None,
        vep_intervals_num: int = None,
    ):
        utils.setup_env()

        # Mocking elastic search password for the full dry run test:
        seqr_loader._read_es_password = Mock(return_value='TEST')

        pipeline = Pipeline(
            namespace=Namespace.TEST,
            name=self._testMethodName,
            analysis_dataset_name=utils.DATASET,
            check_intermediates=False,
            check_expected_outputs=False,
            storage_provider=UnittestStorageProvider(self.out_bucket),
            first_stage=first_stage,
            last_stage=last_stage,
            version=self.timestamp,
            refs=CpgRefData(),
            images=CpgImages(),
            config=dict(
                realignment_shards_num=realignment_shards_num
                or self.realignment_shards_num,
                hc_intervals_num=hc_intervals_num or self.hc_intervals_num,
                jc_intervals_num=jc_intervals_num or self.jc_intervals_num,
                vep_intervals_num=vep_intervals_num or self.vep_intervals_num,
                intervals_path=intervals_path,
            ),
        )
        self.datasets = [pipeline.create_dataset(utils.DATASET)]
        for ds in self.datasets:
            for s_id in self.sample_ids:
                s = ds.add_sample(s_id, s_id)
                s.alignment_input_by_seq_type[utils.SEQ_TYPE] = AlignmentInput(
                    CramPath(utils.TOY_CRAM_BY_SID[s.id]),
                    utils.SEQ_TYPE,
                )
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
            with patch.object(Stage, '_outputs_are_reusable') as mock_reusable:
                mock_reusable.return_value = False
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

        self.assertEqual(
            _cnt('bazam'), len(self.sample_ids) * self.realignment_shards_num
        )
        self.assertEqual(
            _cnt('HaplotypeCaller'), len(self.sample_ids) * self.hc_intervals_num
        )
        self.assertEqual(_cnt('ReblockGVCF'), len(self.sample_ids))
        self.assertEqual(_cnt('GenotypeGVCFs'), self.jc_intervals_num)
        self.assertEqual(_cnt('GenomicsDBImport'), self.jc_intervals_num)
        self.assertEqual(_cnt('MakeSitesOnlyVcf'), self.jc_intervals_num)
        # Indel + SNP create model + SNP scattered
        self.assertEqual(_cnt('VariantRecalibrator'), 2 + self.jc_intervals_num)
        # Twice to each interva: apply indels, apply SNPs
        self.assertEqual(_cnt('ApplyVQSR'), 2 * self.jc_intervals_num)
        # Gather JC, gather siteonly JC, gather VQSR
        self.assertEqual(_cnt('GatherVcfsCloud'), 3)
        self.assertEqual(_cnt('vep '), self.vep_intervals_num)
        self.assertEqual(_cnt('vep_json_to_ht('), 1)
        self.assertEqual(_cnt('annotate_cohort('), 1)
        self.assertEqual(_cnt('annotate_dataset_mt('), len(self.datasets))
        self.assertEqual(_cnt('hailctl dataproc submit'), 1)

    @skip('Running only dry tests in ths module')
    def test_joint_calling(self):
        """
        Stages up to joint calling.
        """
        pipeline = self._setup_pipeline(
            first_stage=GenotypeSample.__name__,
            last_stage=JointGenotyping.__name__,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome1pct/calling_regions.interval_list',
        )
        result = pipeline.run(dry_run=False, wait=True)
        self.assertEqual('success', result.status()['state'])

    @skip('Running only dry tests in ths module. VEP takes too long')
    def test_after_joint_calling(self):
        """
        VQSR needs more inputs than provided by toy regions.
        """
        pipeline = self._setup_pipeline(
            first_stage=Vqsr.__name__,
            last_stage=AnnotateDataset.__name__,
            intervals_path=utils.BASE_BUCKET
            / 'inputs/exome5pct/calling_regions.interval_list',
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
