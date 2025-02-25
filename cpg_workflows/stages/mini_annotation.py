#!/usr/bin/env python3

"""
This file describes a workflow to apply minimal annotations to a callset
Starting from first principles -

- single-sample VCFs
- publicly available reference data
- no specialized tools (or just tools that can be built from source easily)
"""

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.targets import SequencingGroup, MultiCohort, Cohort
from cpg_workflows.utils import ExpectedResultT
from cpg_workflows.workflow import MultiCohortStage, StageInput, StageOutput, stage, CohortStage, SequencingGroupStage, get_multicohort
from cpg_workflows.jobs.bcftools import strip_gvcf_to_vcf, merge_ss_vcfs


@stage
class StripSingleSampleGvcf(SequencingGroupStage):
    """
    Merge all the single-sample VCFs into a single VCF
    do this by finding the gVCFs and stripping them to VCF
    This isn't perfect, it's just to stress test this type of workflow
    We get the same results by condensing the combiner output to joint-VCF
    """

    def expected_outputs(self, sg: SequencingGroup) -> Path:
        return self.tmp_prefix / f'{sg.id}_stripped.vcf.bgz'

    def queue_jobs(self, sg: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        Take the gVCF for this sample and strip it to a VCF (kinda...)
        """

        outputs = self.expected_outputs(sg)

        job = strip_gvcf_to_vcf(
            gvcf=sg.gvcf,
            output=str(outputs),
        )
        return self.make_outputs(sg, data=outputs, jobs=job)


@stage(required_stages=StripSingleSampleGvcf)
class MergeSingleSampleVcfs(CohortStage):
    """
    Merge all the single-sample VCFs into a single VCF
    use a specific region
    """
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return self.tmp_prefix / 'merged_ss_vcfs.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Merge all the single-sample VCFs into a single VCF
        """
        outputs = self.expected_outputs(cohort)
        stripped_vcfs = inputs.as_path_by_target(StripSingleSampleGvcf)

        # get all of those relevant to this cohort
        vcfs_to_merge = [str(stripped_vcfs[sg.id]) for sg in cohort.get_sequencing_groups()]

        job = merge_ss_vcfs(
            input_list=vcfs_to_merge,
            output_file=str(outputs),
            missing_to_ref=True,  # heck, why not
            region_file=config_retrieve(['annotations', 'region']),
        )

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=MergeSingleSampleVcfs)
class AnnotateGnomadFrequenciesWithEchtvar(CohortStage):
    """
    Annotate this cohort with gnomad frequencies
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to generate a matrix table
        """
        return self.tmp_prefix / 'gnomad_frequency_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        apply pre-encoded gnomad frequencies to the merged vcf
        """
        outputs = self.expected_outputs(cohort)

        merged_vcf = get_batch().read_input(str(inputs.as_path(cohort, MergeSingleSampleVcfs)))

        gnomad_annotations = get_batch().read_input(config_retrieve(['annotations', 'echtvar']))

        job = get_batch().new_job('Annotate gnomad frequencies with echtvar')
        job.image(image_path('echtvar'))
        job.command(f'echtvar anno -e {gnomad_annotations} {merged_vcf} {job.output}')
        job.storage('30Gi')
        job.memory('highmem')
        job.cpu(4)

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=[AnnotateGnomadFrequenciesWithEchtvar])
class AnnotateConsequenceWithBcftools(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.tmp_prefix / 'consequence_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        output = self.expected_outputs(cohort)
        gnomad_annotated_vcf = get_batch().read_input(str(inputs.as_path(cohort, AnnotateGnomadFrequenciesWithEchtvar)))

        # get the GFF3 file required to generate consequences
        gff3_file = get_batch().read_input(config_retrieve(['annotations', 'gff3']))

        # get the fasta
        fasta = get_batch().read_input(config_retrieve(['workflow', 'ref_fasta']))

        job = get_batch().new_job('bcftools csq')
        job.image(image_path('bcftools_120'))
        job.cpu(4)
        job.memory('highmem')
        job.storage('10G')

        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        job.command(
            f"""
            bcftools index -t {gnomad_annotated_vcf}
            bcftools csq --force -f {fasta} \
                --local-csq \
                -g {gff3_file} \
                -Oz -o {job.output["vcf.bgz"]} \
                {gnomad_annotated_vcf}
            bcftools index -t {job.output["vcf.bgz"]}
            """
        )

        get_batch().write_output(job.output, output.removesuffix('.vcf.bgz'))

        return self.make_outputs(cohort, data=output, jobs=job)


@stage
class WgetAlphaMissenseTsv(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> ExpectedResultT:
        return to_path(config_retrieve(['storage', 'common', 'analysis'])) / 'reannotation' / 'AlphaMissense38.tsv.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)

        job = get_batch().new_job('wget Alpha missense tsv')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'wget -O {job.output} https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz')
        job.storage('10Gi')

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=WgetAlphaMissenseTsv)
class ReformatAlphaMissenseTsv(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return to_path(config_retrieve(['storage', 'common', 'analysis'])) / 'reannotation' / 'AlphaMissense38.ht.tar.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)

        am_tsv = get_batch().read_input(str(inputs.as_path(multicohort, WgetAlphaMissenseTsv)))

        job = get_batch().new_job('Reformat AlphaMissense TSV')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'convert_alpha_missense --am_tsv {am_tsv} --ht_out AlphaMissense38.ht')
        job.command(f'tar -czf {job.output} AlphaMissense38.ht')
        job.storage('10Gi')
        job.cpu(4)
        job.memory('highmem')
        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[ReformatAlphaMissenseTsv, AnnotateConsequenceWithBcftools])
class CombineAnnotatedVcfAndAlphaMissenseIntoHT(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.prefix / 'annotated_for_reanalysis.ht.tar.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        # convert_vcf_to_mt
        output = self.expected_outputs(cohort)

        # pull the alphamissense TarBall
        alpha_missense_tar = get_batch().read_input(str(inputs.as_path(get_multicohort(), ReformatAlphaMissenseTsv)))

        # get the annotated VCF & index
        vcf = str(inputs.as_path(cohort, AnnotateConsequenceWithBcftools))
        vcf_in = get_batch().read_input_group(**{'vcf.bgz': vcf, 'vcf.bgz.tbi': f'{vcf}.tbi'})['vcf.bgz']

        # bed with region and gene mapping
        region_file = get_batch().read_input(config_retrieve(['annotations', 'region']))

        job = get_batch().new_job('Combine annotated VCF and AlphaMissense data')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'tar -xf {alpha_missense_tar}')
        job.cpu(4)
        job.memory('highmem')
        job.command(
            f'convert_vcf_to_mt '
            f'--input {vcf_in} '
            f'--am alpha_missense.ht '
            f'--gene_bed {region_file} '
            f'--output tmp.mt'
        )
        job.command(f'tar -czf {job.output} tmp.mt')

        # write the output
        get_batch().write_output(job.output, str(output))

        return self.make_outputs(cohort, data=output, jobs=job)

