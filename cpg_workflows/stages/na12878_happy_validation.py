#!/usr/bin/env python3

"""
Content relating to the hap.py validation process, specific to NA12878 from gVCF
"""

from cpg_utils import Path, config, hail_batch

from cpg_utils.config import config_retrieve, dataset_path
from cpg_workflows.jobs.NA12878_validation import run_happy_on_gvcf
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import SequencingGroupStage, StageInput, StageOutput, stage


@stage()
class SingleSampleHappyValidationNA12878(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: SequencingGroup) -> Path:

        return sequencing_group.dataset.prefix() / self.name / f'{sequencing_group.id}.summary.csv'

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:

        # get the input vcf for this sequence group
        input_gvcf = sequencing_group.gvcf.path

        output = self.expected_outputs(sequencing_group)

        job = run_happy_on_gvcf(
            vcf_path=input_gvcf,
            sample_ext_id=sequencing_group.external_id,
            output=output,
        )

        return self.make_outputs(sequencing_group, data=output, jobs=job)
