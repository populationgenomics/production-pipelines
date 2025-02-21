"""
all processes related to echtvar - generating a VCF annotation resource from raw inputs, or appylying those annotations
this is not confined to a specific cohort/project
"""
from os.path import join
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.workflow import stage, MultiCohortStage, StageOutput, StageInput
from cpg_workflows.targets import MultiCohort

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']]
common_folder = join(config_retrieve(['storage', 'common', 'analysis']), 'gnomad', 'echtvar')


@stage
class RunEchtvarOnGnomad(MultiCohortStage):
    """Run echtvar on gnomAD data, generate annotation sources"""
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:

        return {contig: to_path(f'{common_folder}/{contig}.zip') for contig in CANONICAL_CHROMOSOMES}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        run echtvar encode on each gnomadV4 contig
        """
        jobs = []
        outputs = self.expected_outputs(multicohort)
        for contig, output in outputs.items():
            job = get_batch().new_job(f'Run echtvar on {contig}')
            job.image(image_path('echtvar'))

            # some of these are chonky
            job.storage('100G')
            job.cpu(4)
            job.memory('highmem')

            contig_vcf = config_retrieve(['references', 'gnomad_4.1_vcfs', contig])
            contig_in = get_batch().read_input(contig_vcf)
            job.command(f'echtvar encode {job.output} $ECHTVAR_CONFIG {contig_in}')
            get_batch().write_output(job.output, str(output))
            jobs.append(job)

        return self.make_outputs(target=multicohort, jobs=jobs, data=outputs)

