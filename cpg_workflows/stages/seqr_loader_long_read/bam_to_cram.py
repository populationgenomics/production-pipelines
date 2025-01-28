"""
Stage that converts a BAM file to a CRAM file.
Intended for use with long-read BAM files from PacBio.
"""

from cpg_utils import Path, to_path
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
        reference_assembly=config_retrieve(['workflow', 'ref_fasta'], reference_path('broad/ref_fasta')),
    )


@stage(analysis_type=config_retrieve(['workflow', 'bam_to_cram', 'analysis_type'], 'cram'), analysis_keys=['cram'])
class ConvertPacBioBamToCram(SequencingGroupStage):
    """
    Convert a PacBio BAM to a CRAM file.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        if sequencing_group.cram:
            cram_path = sequencing_group.cram
        else:
            cram_path = make_long_read_cram_path(sequencing_group)

        return {'cram': cram_path.path, 'crai': cram_path.index_path}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Using the existing `bam_to_cram` function from the `jobs` module.
        """
        b = get_batch()
        assert (
            sequencing_group.sequencing_type == 'genome'
        ), f'Only genomes are supported for BAM to CRAM conversion, unavailable for {sequencing_group.id}'
        input_bam = b.read_input_group(bam=str(sequencing_group.alignment_input))
        job, output_cram = bam_to_cram.bam_to_cram(
            b=get_batch(),
            input_bam=input_bam,
            extra_label='long_read',
            job_attrs=self.get_job_attrs(sequencing_group),
            requested_nthreads=config_retrieve(['resource_overrides', 'bam_to_cram', 'nthreads'], 1),
            reference_fasta_path=to_path(config_retrieve(['workflow', 'ref_fasta'], reference_path('broad/ref_fasta'))),
            add_rg=config_retrieve(['workflow', 'bam_to_cram', 'add_rg'], False),
        )
        b.write_output(output_cram, str(self.expected_outputs(sequencing_group)['cram']).removesuffix('.cram'))

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=[job])
