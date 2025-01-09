"""
Workflow stages for the ClinvArbitration process

https://github.com/populationgenomics/ClinvArbitration

This takes ClinVar data:
- re-summarises the submitted evidence using new criteria to come to new conclusions
- annotates the resultant VCF representation of the results
- processes the annotated data to create a PM5 resource
    - identify where a missense occurs at the same codon as known pathogenic missense
    - this is used in the Talos pipeline
"""

# mypy: ignore_errors


from datetime import datetime
from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job
from cpg_workflows import get_batch
from cpg_workflows.jobs.simple_vep_annotation import split_and_annotate_vcf
from cpg_workflows.workflow import MultiCohort, MultiCohortStage, StageInput, StageOutput, stage

DATE_STRING: str = datetime.now().strftime('%y-%m')


@stage
class CopyLatestClinvarFiles(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        common_folder = join(config_retrieve(['storage', 'common', 'analysis']), 'clinvarbitration', DATE_STRING)
        return {
            'submission_file': to_path(join(common_folder, 'submission_summary.txt.gz')),
            'variant_file': to_path(join(common_folder, 'variant_summary.txt.gz')),
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        run a wget copy of the relevant files into GCP
        """
        bash_job = get_batch().new_bash_job('wget latest ClinVar raw files')
        bash_job.image(config_retrieve(['workflow', 'driver_image']))

        directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
        sub_file = 'submission_summary.txt.gz'
        var_file = 'variant_summary.txt.gz'

        bash_job.command(f'wget -q {directory}{sub_file} -O {bash_job.subs}')
        bash_job.command(f'wget -q {directory}{var_file} -O {bash_job.vars}')

        outputs = self.expected_outputs(mc)

        get_batch().write_output(bash_job.subs, str(outputs['submission_file']))
        get_batch().write_output(bash_job.vars, str(outputs['variant_file']))

        return self.make_outputs(data=outputs, jobs=bash_job, target=mc)


@stage(required_stages=CopyLatestClinvarFiles, analysis_type='custom', analysis_keys=['clinvar_decisions'])
class GenerateNewClinvarSummary(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a couple of files and a HT as Paths
        """
        common_folder = to_path(
            join(config_retrieve(['storage', 'common', 'analysis']), 'clinvarbitration', DATE_STRING),
        )
        return {
            'clinvar_decisions': common_folder / 'clinvar_decisions.ht',
            'snv_vcf': common_folder / 'pathogenic_snvs.vcf.bgz',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        # relatively RAM intensive, short running task
        clinvarbitrate = get_batch().new_job('Run ClinvArbitration Summary')
        clinvarbitrate.image(image_path('clinvarbitration')).memory('highmem').cpu('2')
        authenticate_cloud_credentials_in_job(clinvarbitrate)

        # declare a resource group, leave the HT path implicit
        clinvarbitrate.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        if sites_to_blacklist := config_retrieve(['workflow', 'site_blacklist'], []):
            blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
            blacklist_string = f' -b {blacklist_sites}'
        else:
            blacklist_string = ''

        var_file = get_batch().read_input(str(inputs.as_path(mc, CopyLatestClinvarFiles, 'variant_file')))
        sub_file = get_batch().read_input(str(inputs.as_path(mc, CopyLatestClinvarFiles, 'submission_file')))

        clinvarbitrate.command(
            f'resummary -v {var_file} -s {sub_file} -o {clinvarbitrate.output} --minimal {blacklist_string}',
        )
        clinvarbitrate.command(f'gcloud storage cp -r {clinvarbitrate.output}.ht {str(outputs["clinvar_decisions"])}')

        # selectively copy back some outputs
        get_batch().write_output(clinvarbitrate.output, str(outputs['snv_vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(target=mc, data=outputs, jobs=clinvarbitrate)


@stage(required_stages=GenerateNewClinvarSummary)
class AnnotateClinvarDecisions(MultiCohortStage):
    """
    take the vcf output from the clinvar stage, and apply VEP annotations
    """

    def expected_outputs(self, mc: MultiCohort) -> Path:
        return (
            to_path(join(config_retrieve(['storage', 'common', 'analysis']), 'clinvarbitration', DATE_STRING))
            / 'annotated_snv.vcf.bgz'
        )

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(mc)

        # delegate the splitting, annotation, and re-merging to this existing method
        _out_file, jobs = split_and_annotate_vcf(
            vcf_in=str(inputs.as_path(mc, GenerateNewClinvarSummary, 'snv_vcf')),
            out_vcf=str(outputs),
        )

        return self.make_outputs(target=mc, jobs=jobs, data=outputs)


@stage(required_stages=AnnotateClinvarDecisions)
class PM5TableGeneration(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a single HT
        """
        common_folder = to_path(
            join(config_retrieve(['storage', 'common', 'analysis']), 'clinvarbitration', DATE_STRING),
        )
        return {
            'clinvar_pm5': common_folder / 'clinvar_pm5.ht',
            'pm5_json': common_folder / 'clinvar_pm5.json',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:

        # declare a resouce group, but don't copy the whole thing back
        clinvarbitrate_pm5 = get_batch().new_job('Run ClinvArbitration PM5')
        clinvarbitrate_pm5.image(image_path('clinvarbitration'))
        authenticate_cloud_credentials_in_job(clinvarbitrate_pm5)

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        vcf = str(inputs.as_path(mc, AnnotateClinvarDecisions))
        annotated_vcf = get_batch().read_input_group(**{'vcf.gz': vcf, 'vcf.gz.tbi': vcf + '.tbi'})['vcf.gz']

        # using a declared resource group and only exporting part of it failed... not sure why
        clinvarbitrate_pm5.command(f'pm5_table -i {annotated_vcf} -o output')
        clinvarbitrate_pm5.command(f'mv output.json {clinvarbitrate_pm5.output}')

        # recursive copy of the HT
        clinvarbitrate_pm5.command(f'gcloud storage cp -r output.ht {str(outputs["clinvar_pm5"])}')

        # also copy back the JSON file
        get_batch().write_output(clinvarbitrate_pm5.output, str(outputs['pm5_json']))

        return self.make_outputs(target=mc, data=outputs, jobs=clinvarbitrate_pm5)


@stage(
    required_stages=[GenerateNewClinvarSummary, AnnotateClinvarDecisions, PM5TableGeneration],
    analysis_type='custom',
)
class PackageForRelease(MultiCohortStage):
    """
    takes the data created so far, and packages it up for release
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:

        common_folder = to_path(
            join(
                config_retrieve(['storage', 'common', 'analysis']),
                'clinvarbitration',
                DATE_STRING,
            ),
        )
        return common_folder / 'clinvarbitration.tar.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Localise all the previously generated data into a folder
        tarball it, and write out as a single file
        """
        tar_output = self.expected_outputs(multicohort)

        # find paths to the previous outputs
        vcf = inputs.as_path(multicohort, AnnotateClinvarDecisions)
        pm5 = inputs.as_dict(multicohort, PM5TableGeneration)
        clinvar_decisions = inputs.as_path(multicohort, GenerateNewClinvarSummary, 'clinvar_decisions')

        job = get_batch().new_job('package clinvarbitration for release')
        job.storage('10G')
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.command(
            f"""
            mkdir clinvarbitration_data
            gcloud storage cp -r {str(pm5["clinvar_pm5"])} clinvarbitration_data
            gcloud storage cp {str(pm5["pm5_json"])} clinvarbitration_data
            gcloud storage cp -r {str(clinvar_decisions)} clinvarbitration_data
            gcloud storage cp "{vcf}*" clinvarbitration_data
            tar -czf {job.output} clinvarbitration_data
        """,
        )
        get_batch().write_output(job.output, str(tar_output))
        return self.make_outputs(multicohort, data=tar_output, jobs=job)
