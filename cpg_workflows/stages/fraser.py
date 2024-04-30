"""
Perform aberrant splicing analysis with FRASER.
"""

from cpg_utils import Path
from cpg_workflows import get_batch
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.jobs import fraser
from cpg_workflows.stages.trim_align import TrimAlignRNA
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    StageInput,
    StageOutput,
    stage,
)


@stage(
    required_stages=TrimAlignRNA,
)
class Fraser(CohortStage):
    """
    Perform aberrant splicing analysis with FRASER.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Generate FRASER outputs.
        """
        dataset_prefix = cohort.get_sequencing_groups()[0].dataset.prefix()
        return {cohort.name: dataset_prefix / 'fraser' / f'{cohort.name}.fds.tar.gz'}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run FRASER.
        """
        sequencing_groups = cohort.get_sequencing_groups()

        bam_or_cram_inputs: list[tuple[BamPath, None] | tuple[CramPath, Path]] = []
        for sequencing_group in sequencing_groups:
            cram_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'cram')
            crai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'crai')
            input_bam_or_cram: BamPath | CramPath | None = None
            try:
                bam_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bam')
                bai_path = inputs.as_path(sequencing_group, TrimAlignRNA, 'bai')
                input_bam_or_cram = BamPath(bam_path, bai_path)
            except KeyError:
                potential_bam_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam'
                potential_bai_path = sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam.bai'
                if potential_bam_path.exists() and potential_bai_path.exists():
                    input_bam_or_cram = BamPath(potential_bam_path, potential_bai_path)
                else:
                    input_bam_or_cram = CramPath(cram_path, crai_path)
            if isinstance(input_bam_or_cram, (BamPath)):
                bam_or_cram_inputs.append((input_bam_or_cram, None))
            elif isinstance(input_bam_or_cram, (CramPath)):
                bam_or_cram_inputs.append((input_bam_or_cram, potential_bam_path))

        j = fraser.fraser(
            b=get_batch(),
            input_bams_or_crams=bam_or_cram_inputs,
            output_fds_path=list(self.expected_outputs(cohort).values())[0],
            cohort_name=cohort.name,
            job_attrs=self.get_job_attrs(),
            overwrite=cohort.forced,
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)
