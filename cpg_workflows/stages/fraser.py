"""
Perform aberrant splicing analysis with FRASER.
"""

import logging
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job
from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows import get_batch
from cpg_workflows.workflow import (
    stage,
    StageInput,
    StageOutput,
    SequencingGroup,
    SequencingGroupStage,
    Cohort,
    CohortStage,
)
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.stages.trim_align import TrimAlignRNA
from cpg_workflows.jobs import bam_to_cram
from cpg_workflows.jobs import fraser


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
        return {
            cohort.name: dataset_prefix / 'fraser' / f'{cohort.name}.output.tar.gz'
        }
    
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Queue a job to run FRASER.
        """
        sequencing_groups = cohort.get_sequencing_groups()
        alignment_inputs: dict[str, dict[str, Path]] = {
            sequencing_group.id: {
                'cram': inputs.as_path(sequencing_group, TrimAlignRNA, 'cram'),
                'crai': inputs.as_path(sequencing_group, TrimAlignRNA, 'crai'),
                'bam': sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam',
                'bai': sequencing_group.dataset.tmp_prefix() / 'bam' / f'{sequencing_group.id}.bam.bai',
            }
            for sequencing_group in sequencing_groups
        }
        cram_bam_inputs: dict[str, dict[str, (BamPath | CramPath)]] = {
            id: {
                'bam': BamPath(aln['bam'], aln['bai']),
                'cram': CramPath(aln['cram'], aln['crai']),
            }
            for id, aln in alignment_inputs.items()
        }
        cram_to_bam_inputs: dict[str, dict[str, (BamPath | CramPath)]] = {
            id: inpt
            for id, inpt in cram_bam_inputs.items()
            if not (alignment_inputs[id]['bam'].exists() and alignment_inputs[id]['bai'].exists())
        }
        bam_inputs: list[BamPath | ResourceGroup] = [
            inpt['bam']
            for id, inpt in cram_bam_inputs.items()
            if (
                (alignment_inputs[id]['bam'].exists() and alignment_inputs[id]['bai'].exists()) and
                isinstance(inpt['bam'], BamPath)
            )
        ]
        jobs: list[Job] = []
        for inpt in cram_to_bam_inputs.values():
            cram = inpt['cram']
            bam = inpt['bam']
            assert isinstance(cram, CramPath)
            assert isinstance(bam, BamPath)
            logging.info(f'Converting {cram} to BAM')
            j, output_bam = bam_to_cram.cram_to_bam(
                b=get_batch(),
                input_cram=cram,
                output_bam=bam,
                job_attrs=self.get_job_attrs(),
                overwrite=cohort.forced,
            )
            jobs.append(j)
            assert isinstance(output_bam, ResourceGroup)
            bam_inputs.append(output_bam)
        j = fraser.fraser(
            b=get_batch(),
            input_bams=bam_inputs,
            output_path=list(self.expected_outputs(cohort).values())[0],
            cohort_name=cohort.name,
            job_attrs=self.get_job_attrs(),
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=j)
