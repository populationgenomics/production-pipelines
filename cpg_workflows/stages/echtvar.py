"""
all processes related to echtvar - generating a VCF annotation resource from raw inputs, or applying those annotations
this is not confined to a specific cohort/project
"""

from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']]
common_folder = join(config_retrieve(['storage', 'common', 'analysis']), 'gnomad', 'echtvar')


@stage
class RunEchtvarOnGnomad(MultiCohortStage):
    """Run echtvar on gnomAD data, generate annotation sources"""

    def expected_outputs(self, multicohort: MultiCohort) -> Path:

        return to_path(f'{common_folder}/gnomad_4.1_echtvar.zip')

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        run echtvar encode on all gnomadV4 contigs combined
        """
        output = self.expected_outputs(multicohort)

        contig_files = []
        for contig in CANONICAL_CHROMOSOMES:
            contig_files.append(get_batch().read_input(config_retrieve(['references', 'gnomad_4.1_vcfs', contig])))

        job = get_batch().new_job('Run echtvar on gnomad v4.1')
        job.image(image_path('echtvar'))

        # some of these are chonky
        job.storage('1000G')
        job.cpu(4)
        job.memory('highmem')

        job.command(f'echtvar encode {job.output} $ECHTVAR_CONFIG {" ".join(contig_files)}')
        get_batch().write_output(job.output, str(output))
        return self.make_outputs(target=multicohort, jobs=job, data=output)
