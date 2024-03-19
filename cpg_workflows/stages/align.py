"""
Stage that generates a CRAM file.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs import align
from cpg_workflows.jobs.align import MissingAlignmentInputException
from cpg_workflows.workflow import (
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(analysis_type='cram', analysis_keys=['cram'])
class Align(SequencingGroupStage):
    """
    Align or re-align input data to produce a CRAM file
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        if sequencing_group.cram:
            cram_path = sequencing_group.cram
        else:
            cram_path = sequencing_group.make_cram_path()

        return {
            'cram': cram_path.path,
            'crai': cram_path.index_path,
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
        """
        # We are not listing markduplicates metrics in expected_outputs:
        # if CRAM exist but metrics are for some reason missing, we want to avoid
        # triggering re-running the expensive alignment again. Samtools stats cover
        # the deduplication rate metric. However, we still want to save markduplicates
        # metrics just in case, as they are more detailed.
        markdup_metrics_path = (
            sequencing_group.dataset.prefix()
            / 'qc'
            / 'markduplicates_metrics'
            / f'{sequencing_group.id}.markduplicates-metrics'
        )

        try:
            jobs = align.align(
                b=get_batch(),
                sequencing_group=sequencing_group,
                output_path=sequencing_group.make_cram_path(),
                job_attrs=self.get_job_attrs(sequencing_group),
                overwrite=sequencing_group.forced,
                out_markdup_metrics_path=markdup_metrics_path,
            )
            return self.make_outputs(
                sequencing_group,
                data=self.expected_outputs(sequencing_group),
                jobs=jobs,
            )
        except MissingAlignmentInputException:
            if get_config()['workflow'].get('skip_sgs_with_missing_input'):
                logging.error(f'No alignment inputs, skipping sequencing group {sequencing_group}')
                sequencing_group.active = False
                return self.make_outputs(sequencing_group, skipped=True)  # return empty output
            else:
                return self.make_outputs(target=sequencing_group, error_msg='No alignment input found')
