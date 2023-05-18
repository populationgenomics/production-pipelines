import unittest
from cpg_workflows.workflow import Workflow, get_workflow
from cpg_workflows.stages.align import Align
from cpg_workflows.stages.genotype import Genotype
from cpg_workflows.stages.joint_genotyping import JointGenotyping
from cpg_workflows.stages.cram_qc import CramMultiQC
from cpg_workflows.stages.gvcf_qc import GvcfMultiQC
from cpg_workflows.stages.large_cohort import LoadVqsr, Frequencies, AncestryPlots
from cpg_utils.config import set_config_paths, get_config

import tempfile
import shutil
import os

# TODO: Remove the hardcoding here

TOML = f"""
[workflow]
dataset_gcp_project = 'fewgenomes'
access_level = 'test'
dataset = 'fewgenomes'
sequencing_type = 'genome'
skip_samples = ['CPG280164']

check_inputs = false
check_intermediates = false
check_expected_outputs = false
path_scheme = 'local'

[storage.default]
analysis = "gs://cpg-fewgenomes-test-analysis"
default = "gs://cpg-fewgenomes-test"
tmp = "gs://cpg-fewgenomes-test-tmp"
upload = "gs://cpg-fewgenomes-test-upload"
web = "gs://cpg-fewgenomes-test-web"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[storage.fewgenomes]
analysis = "gs://cpg-fewgenomes-test-analysis"
default = "gs://cpg-fewgenomes-test"
tmp = "gs://cpg-fewgenomes-test-tmp"
upload = "gs://cpg-fewgenomes-test-upload"
web = "gs://cpg-fewgenomes-test-web"
web_url = "https://test-web.populationgenomics.org.au/fewgenomes"

[hail]
billing_project = 'fewgenomes'
delete_scratch_on_exit = false
backend = 'local'

"""


class TestStageBuild(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.environ["GOOGLE_CLOUD_PROJECT"] = "fewgenomes"

        with tempfile.NamedTemporaryFile(
            mode='w', delete=False, suffix='.toml'
        ) as toml_file:
            toml_file.write(TOML)
            toml_filename = toml_file.name

        config_paths = [toml_filename, 'configs/dev.toml']

        workflow = 'large_cohort'
        wfl_conf_path = f'configs/defaults/{workflow}.toml'
        config_paths = os.environ['CPG_CONFIG_PATH'].split(',') + list(config_paths)
        set_config_paths(config_paths[:1] + [str(wfl_conf_path)] + config_paths[1:])

        self.workflow = get_workflow(dry_run=True)

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_all_stages(self):
        """
        Test regular stage building
        """
        STAGES = [Align, Genotype, CramMultiQC, GvcfMultiQC, JointGenotyping]

        self.workflow.run(stages=STAGES)
        all_stages = []
        for stage in self.workflow.queued_stages:
            all_stages.extend(stage.name for stage in stage.required_stages)

        all_stages = set(
            all_stages + [stage.name for stage in self.workflow.queued_stages]
        )
        self.assertEqual(len(self.workflow.queued_stages), len(all_stages))

    def test_more_stages(self):
        """
        Test regular stage building
        """

        STAGES = [LoadVqsr, Frequencies, AncestryPlots, GvcfMultiQC, CramMultiQC]

        self.workflow.run(stages=STAGES)
        all_stages = []
        for stage in self.workflow.queued_stages:
            all_stages.extend(stage.name for stage in stage.required_stages)

        all_stages = set(
            all_stages + [stage.name for stage in self.workflow.queued_stages]
        )
        print(self.workflow.queued_stages)
        print(all_stages)
        self.assertEqual(len(self.workflow.queued_stages), len(all_stages))

    def test_first_stage(self):
        """
        Test setting the first_stage paramater
        """

        STAGES = [Align, Genotype, CramMultiQC, GvcfMultiQC, JointGenotyping]

        NEW_TOML = f"""
        [workflow]
        first_stages = ['Genotype']

        """
        with tempfile.NamedTemporaryFile(
            mode='w', delete=False, suffix='.toml'
        ) as new_toml_file:
            new_toml_file.write(NEW_TOML)
            new_toml_filename = new_toml_file.name

        with tempfile.NamedTemporaryFile(
            mode='w', delete=False, suffix='.toml'
        ) as toml_file:
            toml_file.write(TOML)
            toml_filename = toml_file.name

        config_paths = [toml_filename, 'configs/dev.toml', new_toml_filename]

        workflow = 'large_cohort'

        config_paths = os.environ['CPG_CONFIG_PATH'].split(',') + list(config_paths)
        self.workflow = get_workflow(dry_run=True)

        set_config_paths(config_paths=config_paths)

        self.workflow.run(stages=STAGES)

        all_stages = []
        for stage in self.workflow.queued_stages:
            all_stages.extend(stage.name for stage in stage.required_stages)

        all_stages = set(
            all_stages + [stage.name for stage in self.workflow.queued_stages]
        )

        all_stages.remove('Align')
        self.assertEqual(len(self.workflow.queued_stages), len(all_stages))
