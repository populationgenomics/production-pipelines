"""
all processes related to echtvar - generating a VCF annotation resource from raw inputs, or applying those annotations
this is not confined to a specific cohort/project
"""

from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import MultiCohort
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, exists, stage

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'whole_genome']]
common_folder = join(config_retrieve(['storage', 'common', 'analysis']), 'gnomad', 'echtvar')
bcftools_image = image_path('bcftools_120')
echtvar_image = image_path('echtvar')


def storage_with_buffer(file_path: str, buffer: int = 10) -> int:
    """
    determine the storage requirement for a file, adding a buffer
    Args:
        file_path ():
        buffer (int): number of GiB to add to the storage requirement
    """
    # determine exact storage requirement, add a buffer for safety and outputs
    return (to_path(file_path).stat().st_size // 1024**3) + buffer


@stage
class RunEchtvarOnGnomad(MultiCohortStage):
    """Run echtvar on gnomAD data, generate annotation sources"""

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        # run for each chromosome, and again for the whole genome
        return {chrom: to_path(f'{common_folder}/gnomad_4.1_{chrom}_trimmed.zip') for chrom in CANONICAL_CHROMOSOMES}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        run echtvar encode on all gnomadV4 contigs, separately and combined
        I'm popping in a quick region filter, to reduce the size of the output files and processing time
        we need to do this once ever, estimated cost $5
        """

        output = self.expected_outputs(multicohort)

        # read this in, once
        bed_file = get_batch().read_input(config_retrieve(['ensembl', 'bed_file']))

        jobs = []
        contig_files = []
        for contig in CANONICAL_CHROMOSOMES:
            # don't do this for the whole genome output
            if contig == 'whole_genome':
                continue

            if exists(output[contig]):
                get_logger().info(f'Skipping echtvar on {contig}, output already exists')
                continue

            # localise this one file
            contig_path = config_retrieve(['references', 'gnomad_4.1_vcfs', contig])
            job_storage = storage_with_buffer(contig_path)
            contig_localised = get_batch().read_input(contig_path)

            trim_job = get_batch().new_bash_job(f'Trim {contig} to Ensembl regions')
            trim_job.storage(f'{job_storage}Gi')
            trim_job.cpu(4)
            trim_job.command(f'bcftools view -R {bed_file} --regions-overlap 2 {contig_localised} -O z -o {trim_job.output}')

            # use the trimmed version for the whole genome
            contig_files.append(trim_job.output)

            # encode that
            contig_job = get_batch().new_job(f'Run echtvar on gnomad v4.1, {contig}')
            contig_job.image(image_path('echtvar'))
            contig_job.storage(f'{job_storage}Gi')
            contig_job.cpu(4)
            contig_job.memory('highmem')

            # run the echtvar encode command
            contig_job.command(f'echtvar encode {contig_job.output} $ECHTVAR_CONFIG {trim_job.output}')
            get_batch().write_output(contig_job.output, str(output[contig]))
            jobs.append(contig_job)

        if not exists(output['whole_genome']):
            job = get_batch().new_job('Run echtvar on gnomad v4.1, whole genome')
            job.image(image_path('echtvar'))
            job.storage('100G')
            job.cpu(4)
            job.memory('highmem')
            job.command(f'echtvar encode {job.output} $ECHTVAR_CONFIG {" ".join(contig_files)}')
            get_batch().write_output(job.output, str(output['whole_genome']))
            jobs.append(job)

        return self.make_outputs(target=multicohort, jobs=jobs, data=output)
