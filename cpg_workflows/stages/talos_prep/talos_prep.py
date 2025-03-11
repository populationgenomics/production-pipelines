"""
creates a new workflow designed to take raw calls and apply quick annotation
this workflow is designed to be something which could be executed easily off-site by non-CPG users
"""
from cpg_utils import to_path, Path
from cpg_workflows.workflow import get_workflow, get_multicohort, get_batch, Stage, stage, MultiCohortStage, \
    CohortStage, StageInput, StageOutput
from cpg_workflows.targets import MultiCohort, Dataset, Cohort

from cpg_workflows.stages.talos import query_for_latest_hail_object


@stage
class QuickAnnotation(CohortStage):
    """
    extract some plain calls from a joint-callset
    these calls are a region-filtered subset, limited to genic regions
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        # get prefix for this cohort
        cohort_prefix = get_workflow().cohort_prefix(cohort)
        return {
            'vcf': cohort_prefix / 'extracted.vcf.bgz',
            'sites_only': cohort_prefix / 'extracted.sites_only.vcf.bgz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """

        Args:
            cohort ():
            inputs ():

        Returns:

        """
        ...
