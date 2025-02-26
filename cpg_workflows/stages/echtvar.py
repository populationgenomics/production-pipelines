"""
all processes related to echtvar - generating a VCF annotation resource from raw inputs, or applying those annotations
this is not confined to a specific cohort/project
"""

from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import MultiCohort
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, exists, stage
from cpg_workflows.utils import get_logger

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'whole_genome']]
common_folder = join(config_retrieve(['storage', 'common', 'analysis']), 'gnomad', 'echtvar')


@stage
class RunEchtvarOnGnomad(MultiCohortStage):
    """Run echtvar on gnomAD data, generate annotation sources"""

    def expected_outputs(self, multicohort: MultiCohort) -> Path:

        # run for each chromosome, and again for the whole genome
        return {chrom: to_path(f'{common_folder}/gnomad_4.1_{chrom}.zip') for chrom in CANONICAL_CHROMOSOMES}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        run echtvar encode on all gnomadV4 contigs, separately and combined
        we need to do this once ever, estimated cost $5
        """

        output = self.expected_outputs(multicohort)

        contig_files = []
        for contig in CANONICAL_CHROMOSOMES:
            # don't do this for the whole genome output
            if contig == 'whole_genome':
                continue

            # localise this one file
            contig_localised = get_batch().read_input(config_retrieve(['references', 'gnomad_4.1_vcfs', contig]))
            # add to the list of inputs for the whole genome job
            contig_files.append(contig_localised)

            if exists(output[contig]):
                get_logger().info(f'Skipping echtvar on {contig}, output already exists')
                continue
            # create and resource a job
            contig_job = get_batch().new_job(f'Run echtvar on gnomad v4.1, {contig}')
            contig_job.image(image_path('echtvar'))
            contig_job.storage('10Gi')
            contig_job.cpu(4)
            contig_job.memory('highmem')
            # run the echtvar encode command
            contig_job.command(f'echtvar encode {contig_job.output} $ECHTVAR_CONFIG {contig_localised}')
            get_batch().write_output(contig_job.output, str(output[contig]))

        job = get_batch().new_job('Run echtvar on gnomad v4.1, whole genome')
        job.image(image_path('echtvar'))
        job.storage('700G')
        job.cpu(4)
        job.memory('highmem')

        job.command(f'echtvar encode {job.output} $ECHTVAR_CONFIG {" ".join(contig_files)}')
        get_batch().write_output(job.output, str(output['whole_genome']))

        return self.make_outputs(target=multicohort, jobs=job, data=output)
