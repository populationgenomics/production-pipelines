"""
This workflow is designed to generate a minimally annotated dataset for use in Talos
this workflow is designed to be something which could be executed easily off-site by non-CPG users

- We generate a VCF from the cohort's joint-call matrixtable
    - other sites would accomplish the same through a VCF merge
    - this VCF is region limited to Ensembl genes +/- 2kb
    - ROI generated here https://github.com/populationgenomics/references/blob/main/reference_generating_scripts/generate_bed_from_ensembl_manager.sh
    - Talos analysis is gene-based, so region-filtering reduces size and increases speed of everything downstream
    - benchmarking will be useful in adjusting this if we need to

- During this extraction from the MT all INFO fields are dropped, to make sure we have a blank slate

- We generate a sites-only version of the same VCF to annotate (less data to pass around)

- We remove any non-PASS sites from the MT prior to export
    - This is necessary because we don't pull the VQSR headers when applying QC annotations, just the filters
    - >> "FILTER 'VQSRTrancheSNP99.90to100.00' is not defined in the header"
    - There's a longer-term fix for this, but it's not important to this workflow

- We annotate the sites-only VCF with gnomad frequencies
    - uses the gnomAD 4.1 joint exomes + genomes dataset

- We annotate the sites-only VCF with consequences using BCFtools

- We reformat the annotations into a HailTable
    - Split up the CSQ string into a Hail Array of Hail Structs
    - integrate MANE identifiers per-transcript
    - integrate AlphaMissense scores per-transcript

- Load the annotations Table and the full VCF - hop annotations into the full MT
"""

from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_workflows.jobs.gcloud_composer import gcloud_compose_vcf_from_manifest
from cpg_workflows.stages.talos import query_for_latest_hail_object
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import get_logger, ExpectedResultT
from cpg_workflows.workflow import CohortStage, StageInput, StageOutput, get_batch, get_workflow, stage


SHARD_MANIFEST = 'shard-manifest.txt'


@stage
class ExtractVcfFromCohortMt(CohortStage):
    """
    extract some plain calls from a joint-callset
    these calls are a region-filtered subset, limited to genic regions
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        # get prefix for this cohort
        cohort_prefix = get_workflow().cohort_prefix(cohort, category='tmp')

        return {
            # this will be the write path for fragments of sites-only VCF
            'vcf_dir': str(cohort_prefix / f'{cohort.id}.vcf.bgz'),
            # this will be the file which contains the name of all fragments
            'vcf_manifest': cohort_prefix / f'{cohort.id}.vcf.bgz' / SHARD_MANIFEST,
            # this will be the write path for fragments of sites-only VCF
            'sites_only_vcf_dir': str(cohort_prefix / f'{cohort.id}_separate.vcf.bgz'),
            # this will be the file which contains the name of all fragments
            'sites_only_vcf_manifest': cohort_prefix / f'{cohort.id}_separate.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        script is called extract_vcf_from_mt
        """

        # either get a mt from metamist, or take one from config
        if (input_mt := config_retrieve(['workflow', cohort.id, 'mt'], None)) is None:
            get_logger().info(f'No config entry retrieved at workflow/{cohort.name}/mt')
            input_mt = query_for_latest_hail_object(
                dataset=cohort.analysis_dataset.name,
                analysis_type='matrixtable',
                object_suffix='.mt',
            )

        outputs = self.expected_outputs(cohort)

        # get the BED file - does not need to be localised
        ensembl_version = config_retrieve(['workflow', 'ensembl_version'], 113)
        bed = reference_path(f'ensembl_{ensembl_version}/merged_bed')

        job = get_batch().new_job(f'Extract VCF representations from {input_mt} for {cohort.id}')
        job.storage('10Gi')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'extract_vcf_from_mt --mt {input_mt} --out {str(outputs)} --bed {bed}')

        return self.make_outputs(cohort, outputs, jobs=job)


@stage(required_stages=ExtractVcfFromCohortMt)
class ConcatenateVcfFragments(CohortStage):
    """
    glue all the mini-VCFs together
    """
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:

        cohort_prefix = get_workflow().cohort_prefix(cohort)
        return {
            'vcf': cohort_prefix / f'{cohort.id}_reassembled.vcf.bgz',
            'sites_only': cohort_prefix / f'{cohort.id}_sites_only_reassembled.vcf.bgz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        trigger a rolling merge using gcloud compose, gluing all the individual files together
        """

        outputs = self.expected_outputs(cohort)
        prior_outputs = inputs.as_dict(cohort, ExtractVcfFromCohortMt)

        full_manifest = prior_outputs['vcf_manifest']
        sites_manifest = prior_outputs['sites_only_vcf_manifest']

        all_jobs = []

        for in_path, out_path in [(full_manifest, outputs['vcf']), (sites_manifest, outputs['sites_only'])]:
            if not in_path.exists():
                raise ValueError(
                    f'Manifest file {str(in_path)} does not exist, '
                    f'run the rd_combiner workflow with workflows.last_stages=[ExtractVcfFromCohortMt]',
                )
            jobs = gcloud_compose_vcf_from_manifest(
                manifest_path=in_path,
                intermediates_path=str(self.tmp_prefix / 'temporary_compose_intermediates'),
                output_path=str(out_path),
                job_attrs={'stage': self.name},
            )

            all_jobs.extend(jobs)
        return self.make_outputs(cohort, data=outputs, jobs=all_jobs)


@stage(required_stages=ConcatenateVcfFragments)
class AnnotateGnomadFrequenciesWithEchtvar(CohortStage):
    """
    Annotate this cohort joint-call VCF with gnomad frequencies, write to tmp storage
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return self.get_stage_cohort_prefix(cohort=cohort, category='tmp') / 'gnomad_frequency_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(cohort)

        sites_vcf = get_batch().read_input(str(inputs.as_dict(cohort, ConcatenateVcfFragments)['sites_only']))

        # this is a single whole-genome file, generated by the echtvar workflow
        gnomad_annotations = get_batch().read_input(config_retrieve(['annotations', 'echtvar']))

        job = get_batch().new_job('Annotate gnomad frequencies with echtvar')
        job.image(image_path('echtvar'))
        job.command(f'echtvar anno -e {gnomad_annotations} {sites_vcf} {job.output}')
        job.storage('20Gi')
        job.memory('highmem')
        job.cpu(4)

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage(required_stages=AnnotateGnomadFrequenciesWithEchtvar)
class AnnotateConsequenceWithBcftools(CohortStage):
    """
    Take the VCF with gnomad frequencies, and annotate with consequences using BCFtools
    Writes into a cohort-specific permanent folder
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort=cohort, category='tmp') / 'consequence_annotated.vcf.bgz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        output = self.expected_outputs(cohort)
        gnomad_annotated_vcf = get_batch().read_input(str(inputs.as_path(cohort, AnnotateGnomadFrequenciesWithEchtvar)))

        # get the GFF3 file required to generate consequences
        ensembl_version = config_retrieve(['workflow', 'ensembl_version'], 113)
        gff3_file = get_batch().read_input(reference_path(f'ensembl_{ensembl_version}/gff3'))

        # get the fasta
        fasta = get_batch().read_input(reference_path('broad/ref_fasta'))

        job = get_batch().new_job('bcftools csq')
        job.image(image_path('bcftools_120'))
        job.cpu(4)
        job.memory('highmem')
        job.storage('20G')

        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # the echtvar image doesn't have a tool to index, so first add that
        # then run the csq command:
        # -g is the GFF3 file
        # -B 10 is where to truncate long protein changes
        # --local-csq indicates we want each variant annotated independently (haplotype unaware)
        job.command(
            f"""
            bcftools index -t {gnomad_annotated_vcf}
            bcftools csq --force -f {fasta} \
                --local-csq \
                -g {gff3_file} \
                -B 10 \
                -Oz -o {job.output["vcf.bgz"]} \
                --write-index=tbi \
                {gnomad_annotated_vcf}
            """,
        )

        get_batch().write_output(job.output, str(output).removesuffix('.vcf.bgz'))

        return self.make_outputs(cohort, data=output, jobs=job)


@stage(required_stages=AnnotateConsequenceWithBcftools)
class ProcessAnnotatedSitesOnlyVcfIntoHt(CohortStage):
    """
    Join the annotated sites-only VCF, with AlphaMissense, and with gene/transcript information
    exporting as a HailTable
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        # output will be a tarball, containing the {cohort.id}_annotations.ht directory
        return self.get_stage_cohort_prefix(cohort=cohort, category='tmp') / f'{cohort.id}_annotations.ht.tar.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        output = self.expected_outputs(cohort)

        # pull the alphamissense TarBall location from config, and localise it
        alphamissense_tar = get_batch().read_input(reference_path('alphamissense/ht_tar'))

        # mane version for gene details
        mane_version = config_retrieve(['workflow', 'mane_version'], '1.4')
        mane_json = get_batch().read_input(reference_path(f'mane_{mane_version}/json'))

        # ensembl version used to generate region of interest
        ensembl_version = config_retrieve(['workflow', 'ensembl_version'], 113)
        gene_roi = get_batch().read_input(reference_path(f'ensembl_{ensembl_version}/bed'))

        # get the annotated VCF & index
        vcf = str(inputs.as_path(cohort, AnnotateConsequenceWithBcftools))
        vcf_in = get_batch().read_input_group(**{'vcf.bgz': vcf, 'vcf.bgz.tbi': f'{vcf}.tbi'})['vcf.bgz']

        job = get_batch().new_job('Combine annotated VCF, AlphaMissense, and MANE annotations')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'tar -xf {alphamissense_tar} -C $BATCH_TMPDIR')
        job.cpu(4)
        job.storage('20Gi')
        job.memory('highmem')
        job.command(
            f'convert_annotated_vcf_to_ht '
            f'--input {vcf_in} '
            f'--am ${{BATCH_TMPDIR}}/alphamissense_38.ht '
            f'--gene_bed {gene_roi} '
            f'--mane {mane_json} '
            f'--output {cohort.id}_annotations.ht',
        )
        job.command(f'tar -czf {job.output} {cohort.id}_annotations.ht')

        # write the output
        get_batch().write_output(job.output, str(output))

        return self.make_outputs(cohort, data=output, jobs=job)


@stage(required_stages=[ProcessAnnotatedSitesOnlyVcfIntoHt, ExtractVcfFromCohortMt])
class JumpAnnotationsFromHtToFinalMt(CohortStage):
    """
    Join the annotated sites-only VCF, with AlphaMissense, and with gene/transcript information
    exporting as a HailTable
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        return self.get_stage_cohort_prefix(cohort=cohort) / f'{cohort.id}_talos_ready.mt.tar.gz'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:

        output = self.expected_outputs(cohort)

        # get the full VCF & index
        vcf = str(inputs.as_dict(cohort, ConcatenateVcfFragments)['vcf'])
        vcf_in = get_batch().read_input_group(**{'vcf.bgz': vcf, 'vcf.bgz.tbi': f'{vcf}.tbi'})['vcf.bgz']

        # get the table of compressed annotations
        annotations = get_batch().read_input(str(inputs.as_path(cohort, ProcessAnnotatedSitesOnlyVcfIntoHt)))

        job = get_batch().new_job('Combine annotations and full VCF')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(f'tar -xf {annotations} -C $BATCH_TMPDIR')
        job.cpu(4)
        job.memory('highmem')
        job.storage('100Gi')
        job.command(
            f'transfer_annotations_to_vcf '
            f'--input {vcf_in} '
            f'--annotations ${{BATCH_TMPDIR}}/{cohort.id}_annotations.ht '
            f'--output {cohort.id}_talos_ready.mt',
        )
        # make some more room on the VM?
        job.command(f'rm -r ${{BATCH_TMPDIR}}/{cohort.id}_annotations.ht')

        # create a compressed tarball
        job.command(f'tar -czf {job.output} {cohort.id}_talos_ready.mt')

        # write the output
        get_batch().write_output(job.output, str(output))

        return self.make_outputs(cohort, data=output, jobs=job)
