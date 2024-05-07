"""
Stage that converts a BAM file to a CRAM file.
"""

from cpg_utils import Path
from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import bam_to_cram

from cpg_workflows.workflow import (
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)

def make_long_read_cram_path(sequencing_group: SequencingGroup) -> CramPath:
    """
    Path to a CRAM file. Not checking its existence here.
    """
    path_prefix = config_retrieve(['workflow', 'long_read_cram_path_prefix'], False)
    if not path_prefix:
        raise ValueError('Missing long_read_path_prefix in the config')
    
    path: Path = sequencing_group.dataset.prefix() / 'cram' / 'long_read' / path_prefix / f'{sequencing_group.id}.cram'
    return CramPath(
        path=path,
        index_path=path.with_suffix('.cram.crai'),
        reference_assembly=reference_path('broad/ref_fasta'),
    )


@stage(analysis_type='cram', analysis_keys=['cram'])
class BamToCram(SequencingGroupStage):
    """
    Convert a BAM to a CRAM file.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        if sequencing_group.cram:
            cram_path = sequencing_group.cram
        else:
            cram_path = make_long_read_cram_path(sequencing_group)

        return {
            'cram': cram_path.path,
            'crai': cram_path.index_path,
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Using the "bam_to_cram" function implemented in the `jobs` module.
        """
        jobs = bam_to_cram.bam_to_cram(
            b=get_batch(),
            input_bam=sequencing_group.alignment_input_by_seq_type.get('genome'),
            extra_label='long_read',
            job_attrs=self.get_job_attrs(sequencing_group),
            requested_nthreads=1,
        )
        return self.make_outputs(
            sequencing_group,
            data=self.expected_outputs(sequencing_group),
            jobs=jobs,
        )
