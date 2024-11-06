"""
Stages that say hello
"""

import json
from functools import lru_cache

from google.api_core.exceptions import PermissionDenied

import hailtop.batch as hb
from hail.utils.java import Env

from cpg_utils import Path, to_path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, image_path, reference_path, try_get_ar_guid
from cpg_utils.hail_batch import get_batch, init_batch
from cpg_workflows.batch import override_jar_spec
from cpg_workflows.jobs import hello
from cpg_workflows.stages.gatk_sv.gatk_sv_common import get_images, get_references, queue_annotate_sv_jobs
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.targets import Cohort, Dataset, MultiCohort, SequencingGroup
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    CohortStage,
    DatasetStage,
    MultiCohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
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

        print('Calling init_batch')
        init_batch(worker_memory='highmem', driver_memory='highmem', driver_cores=4)
        print('Done calling init_batch')

        print('Going to do the thing')
        override_jar_spec(config_retrieve(['workflow', 'jar_spec_revision'], None))
        print('Done the thing')

        bkend = Env.backend()
        print(f'{bkend=}')
        print(bkend.debug_info())

        assert seqgroup.cram
        jobs = hello.hello_job(
            cram_path=seqgroup.cram,
            job_attrs=self.get_job_attrs(seqgroup),
            output_base_path=seqgroup.dataset.prefix() / seqgroup.id,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)
