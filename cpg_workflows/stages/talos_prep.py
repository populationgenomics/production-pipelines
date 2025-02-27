#!/usr/bin/env python3

"""
This file describes a workflow to apply minimal annotations to a callset
Starting from first principles -

- single-sample (g)VCFs
- publicly available reference data
- no specialized tools (or just tools that can be built from source easily)
"""

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.bcftools import merge_ss_vcfs, strip_gvcf_to_vcf
from cpg_workflows.targets import Cohort, MultiCohort, SequencingGroup
from cpg_workflows.workflow import (
    CohortStage,
    MultiCohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    get_multicohort,
    stage,
)

REANNOTATION_DIR = to_path(config_retrieve(['storage', 'common', 'analysis'])) / 'reannotation'


@stage
class WgetEnsemblGffFile(MultiCohortStage):
    """
    Reformat the MANE data into a dictionary format
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        version = config_retrieve(['annotations', 'ensembl_version'])
        return REANNOTATION_DIR / f'ensembl_{version}.gff3.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        version = config_retrieve(['annotations', 'ensembl_version'])
        ensembl_url = config_retrieve(['annotations', 'ensembl_url']).format(v=version)

        job = get_batch().new_job('wget Ensembl GFF3 file')
        job.image(config_retrieve(['workflow', 'driver_image']))

        job.command(f'wget {ensembl_url} -O {job.output}')
        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[WgetEnsemblGffFile])
class GenerateGeneRoi(MultiCohortStage):
    """
    parse the Ensembl GFF3 file, and generate a BED file of gene regions
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        version = config_retrieve(['annotations', 'ensembl_version'])
        return REANNOTATION_DIR / f'ensembl_{version}.bed'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)
        gff3_file = get_batch().read_input(str(inputs.as_path(multicohort, WgetEnsemblGffFile)))

        job = get_batch().new_job('Generate Gene ROI')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'generate_gene_roi --gff3 {gff3_file} --output {job.output}')

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage
class StripSingleSampleGvcf(SequencingGroupStage):
    """
    Merge all the single-sample VCFs into a single VCF
    do this by finding the gVCFs and stripping them to VCF
    This isn't perfect, it's just to stress test this type of workflow
    We get the same results internally by condensing the combiner output to joint-VCF
    """

    def expected_outputs(self, sg: SequencingGroup) -> Path:
        return sg.dataset.prefix(category='tmp') / self.name / f'{sg.id}_stripped.vcf.bgz'

    def queue_jobs(self, sg: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        Take the gVCF for this sample and strip it to a VCF (kinda...)
        """

        outputs = self.expected_outputs(sg)

        if sg.gvcf is None:
            raise ValueError(f'No gVCF for {sg.id}')
        job = strip_gvcf_to_vcf(gvcf=sg.gvcf, output=str(outputs))
        return self.make_outputs(sg, data=outputs, jobs=job)


@stage(required_stages=[StripSingleSampleGvcf, GenerateGeneRoi])
class MergeSingleSampleVcfs(CohortStage):
    """
    Merge all the single-sample VCFs into a single VCF
    use a specific region
    Writes to a cohort-specific temp folder
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort=cohort, category='tmp') / 'merged_ss_vcfs.vcf.bgz'

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
            missing_to_ref=True,
            region_file=str(inputs.as_path(get_multicohort(), GenerateGeneRoi)),
        )

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=MergeSingleSampleVcfs)
class AnnotateGnomadFrequenciesWithEchtvar(CohortStage):
    """
    Annotate this cohort joint-call VCF with gnomad frequencies
    Temp storage, open to reconsideration on this
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return self.get_stage_cohort_prefix(cohort=cohort, category='tmp') / 'gnomad_frequency_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(cohort)

        merged_vcf = get_batch().read_input(str(inputs.as_path(cohort, MergeSingleSampleVcfs)))

        # this is a single whole-genome file, generated by the echtvar workflow
        gnomad_annotations = get_batch().read_input(config_retrieve(['annotations', 'echtvar']))

        job = get_batch().new_job('Annotate gnomad frequencies with echtvar')
        job.image(image_path('echtvar'))
        job.command(f'echtvar anno -e {gnomad_annotations} {merged_vcf} {job.output}')
        job.storage('30Gi')
        job.memory('highmem')
        job.cpu(4)

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=[AnnotateGnomadFrequenciesWithEchtvar, WgetEnsemblGffFile])
class AnnotateConsequenceWithBcftools(CohortStage):
    """
    Take the VCF with gnomad frequencies, and annotate with consequences using BCFtools
    Writes into a cohort-specific permanent folder
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort=cohort) / 'consequence_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        output = self.expected_outputs(cohort)
        gnomad_annotated_vcf = get_batch().read_input(str(inputs.as_path(cohort, AnnotateGnomadFrequenciesWithEchtvar)))

        # get the GFF3 file required to generate consequences
        gff3_file = get_batch().read_input(str(inputs.as_path(get_multicohort(), WgetEnsemblGffFile)))

        # get the fasta
        fasta = get_batch().read_input(config_retrieve(['references', 'broad', 'ref_fasta']))

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
                -B 10 \
                -Oz -o {job.output["vcf.bgz"]} \
                {gnomad_annotated_vcf}
            bcftools index -t {job.output["vcf.bgz"]}
            """,
        )

        get_batch().write_output(job.output, str(output).removesuffix('.vcf.bgz'))

        return self.make_outputs(cohort, data=output, jobs=job)


@stage
class MakeManeJson(MultiCohortStage):
    """
    Download MANE data from Ensembl
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        version = config_retrieve(['annotations', 'mane_version'])
        return {
            'mane_summary': REANNOTATION_DIR / f'mane_{version}.summary.txt.gz',
            'mane_json': REANNOTATION_DIR / f'mane_{version}.json',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)

        mane_version = config_retrieve(['annotations', 'mane_version'])
        mane_url = config_retrieve(['annotations', 'mane_url']).format(v=mane_version)

        job = get_batch().new_job('Get and Reformat MANE summary data')
        job.image(config_retrieve(['workflow', 'driver_image']))

        job.declare_resource_group(output={'summary.txt.gz': '{root}.summary.txt.gz', 'json': '{root}.json'})

        job.command(f'wget {mane_url} -O {job.output["summary.txt.gz"]}')
        job.command(f'reformat_mane_summary --input {job.output["summary.txt.gz"]} --output {job.output["json"]}')
        get_batch().write_output(job.output, str(outputs['mane_json']).removesuffix('.json'))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage
class WgetAlphaMissenseTsv(MultiCohortStage):
    """
    Pull the AlphaMissense TSV file from Zenodo
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return REANNOTATION_DIR / 'AlphaMissense.tsv.gz'

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
    """
    Reformat the AlphaMissense TSV file into a Hail Table
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return REANNOTATION_DIR / 'AlphaMissense.ht.tar.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(multicohort)

        am_tsv = get_batch().read_input(str(inputs.as_path(multicohort, WgetAlphaMissenseTsv)))

        job = get_batch().new_job('Reformat AlphaMissense TSV')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'alphamissense_to_ht --am_tsv {am_tsv} --ht_out AlphaMissense.ht')
        job.command(f'tar -czf {job.output} AlphaMissense.ht')
        job.storage('10Gi')
        job.cpu(4)
        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(
    required_stages=[
        ReformatAlphaMissenseTsv,
        AnnotateConsequenceWithBcftools,
        GenerateGeneRoi,
        MakeManeJson,
    ],
    analysis_type='talos_prep',
)
class CombineAnnotatedVcfAndAlphaMissenseIntoMt(CohortStage):
    """
    Join the annotated VCF, with AlphaMissense, and with gene/transcript information
    exporting as a Hail MatrixTable
    definitely permanent storage for this one
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort=cohort) / 'annotated_for_reanalysis.ht.tar.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        # convert_vcf_to_mt
        output = self.expected_outputs(cohort)

        # pull the alphamissense TarBall
        alphamissense_tar = get_batch().read_input(str(inputs.as_path(get_multicohort(), ReformatAlphaMissenseTsv)))

        mane_json = get_batch().read_input(str(inputs.as_path(get_multicohort(), MakeManeJson, 'mane_json')))

        gene_roi = get_batch().read_input(str(inputs.as_path(get_multicohort(), GenerateGeneRoi)))

        # get the annotated VCF & index
        vcf = str(inputs.as_path(cohort, AnnotateConsequenceWithBcftools))
        vcf_in = get_batch().read_input_group(**{'vcf.bgz': vcf, 'vcf.bgz.tbi': f'{vcf}.tbi'})['vcf.bgz']

        job = get_batch().new_job('Combine annotated VCF and AlphaMissense data')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'tar -xf {alphamissense_tar}')
        job.cpu(4)
        job.memory('highmem')
        job.command(
            f'convert_vcf_to_mt '
            f'--input {vcf_in} '
            f'--am AlphaMissense.ht '
            f'--gene_bed {gene_roi} '
            f'--mane {mane_json} '
            f'--output annotated.mt',
        )
        job.command(f'tar -czf {job.output} annotated.mt')

        # write the output
        get_batch().write_output(job.output, str(output))

        return self.make_outputs(cohort, data=output, jobs=job)
