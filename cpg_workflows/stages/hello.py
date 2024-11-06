"""
Stages that say hello
"""

import hailtop.batch as hb
from hail.utils.java import Env

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.jobs import hello
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.workflow import (
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@stage
class HelloStage(SequencingGroupStage):
    """
    Hello
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'hello': seqgroup.dataset.prefix() / f'{seqgroup.id}.hello',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        #print('Calling init_batch')
        init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)
        #print('Done calling init_batch')

        #print('Going to do the thing')
        if jar_spec := config_retrieve(['workflow', 'jar_spec_revision'], False):
            override_jar_spec(jar_spec)
        #print('Done the thing')

        bkend = Env.backend()
        #print(f'{bkend=}')
        #print(bkend.debug_info())

        assert seqgroup.cram
        jobs = hello.hello_job(
            cram_path=seqgroup.cram,
            job_attrs=self.get_job_attrs(seqgroup),
            output_base_path=seqgroup.dataset.prefix() / seqgroup.id,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)
