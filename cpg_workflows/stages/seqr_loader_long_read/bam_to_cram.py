"""
Stage that converts a BAM file to a CRAM file.
Intended for use with long-read BAM files from PacBio.
"""

from cpg_utils import Path
from cpg_utils.config import config_retrieve, reference_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import bam_to_cram
from cpg_workflows.workflow import SequencingGroup, SequencingGroupStage, StageInput, StageOutput, stage


def make_long_read_cram_path(sequencing_group: SequencingGroup) -> CramPath:
    """
    Path to a CRAM file. Not checking its existence here.
    """
    path: Path = sequencing_group.dataset.prefix() / 'pacbio' / 'cram' / f'{sequencing_group.id}.cram'
    return CramPath(
        path=path,
        index_path=path.with_suffix('.cram.crai'),
        reference_assembly=reference_path('broad/ref_fasta'),
    )


@stage(analysis_type=config_retrieve(['workflow', 'bam_to_cram_analysis_type'], 'cram'), analysis_keys=['cram'])
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

        return {'cram': cram_path.path, 'cram.crai': cram_path.index_path}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Using the existing `bam_to_cram` function from the `jobs` module.
        """
        b = get_batch()
        assert (
            sequencing_group.sequencing_type == 'genome'
        ), f'Only genome sequencing type is supported for BAM to CRAM conversion, unavailable for {sequencing_group.id}'

        outputs = self.expected_outputs(sequencing_group)
        input_bam = b.read_input_group(bam=str(sequencing_group.alignment_input))
        job, output_cram = bam_to_cram.bam_to_cram(
            b=get_batch(),
            input_bam=input_bam,
            extra_label='long_read',
            job_attrs=self.get_job_attrs(sequencing_group),
            requested_nthreads=1,
            reference_fasta_path=reference_path('broad/ref_fasta'),
        )
        b.write_output(output_cram, str(outputs['cram']).removesuffix('.cram'))

        return self.make_outputs(sequencing_group, data=outputs, jobs=[job])
