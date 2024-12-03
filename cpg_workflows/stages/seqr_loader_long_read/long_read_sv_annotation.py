"""
Workflow for finding long-read SV files, updating file contents, merging, and annotating the results
"""

from functools import cache

from google.api_core.exceptions import PermissionDenied

from cpg_utils import Path
from cpg_utils.config import AR_GUID_NAME, config_retrieve, image_path, try_get_ar_guid
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.seqr_loader_sv import annotate_cohort_jobs_sv, annotate_dataset_jobs_sv
from cpg_workflows.stages.gatk_sv.gatk_sv_common import queue_annotate_strvctvre_job, queue_annotate_sv_jobs
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.targets import Dataset, MultiCohort, SequencingGroup
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import (
    DatasetStage,
    MultiCohortStage,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    get_multicohort,
    get_workflow,
    stage,
)
from metamist.graphql import gql, query

VCF_QUERY = gql(
    """
    query MyQuery($dataset: String!) {
  project(name: $dataset) {
    sequencingGroups {
      id
      analyses(type: {eq: "pacbio_vcf"}) {
        output
      }
    }
  }
}""",
)


@cache
def query_for_sv_vcfs(dataset_name: str) -> dict[str, str]:
    """
    query metamist for the PacBio SV VCFs
    return a dictionary of each CPG ID and its corresponding VCF
    this is cached - used in a SequencingGroupStage, but we only want to query for it once instead of once/SG

    Args:
        dataset_name (str):

    Returns:
        a dictionary of the SG IDs and their phased SV VCF
    """
    return_dict: dict[str, str] = {}
    analysis_results = query(VCF_QUERY, variables={'dataset': dataset_name})
    for sg_id_section in analysis_results['project']['sequencingGroups']:
        for analysis in sg_id_section['analyses']:
            if analysis['output'].endswith('.SVs.phased.vcf.gz'):
                return_dict[sg_id_section['id']] = analysis['output']

    return return_dict


@stage
class ReFormatPacBioSVs(SequencingGroupStage):
    """
    take each of the long-read SV VCFs, and re-format the contents
    the huge REF strings are just not required
    a symbolic ALT allele (instead of "N") is required for annotation
    adds a unique ID to each record, so that GATK-SV can do its weird sorting
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.prefix() / 'pacbio' / 'modified_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_reformatted_renamed_lr_svs.vcf.bgz',
            'index': sgid_prefix / f'{sequencing_group.id}_reformatted_renamed_lr_svs.vcf.bgz.tbi',
        }

    def queue_jobs(self, sg: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        a python job to change the VCF contents
        - use a python job to modify all VCF lines, in-line
        - use a follow-up job to block-gzip and index
        """

        # find the vcf for this SG
        query_result = query_for_sv_vcfs(dataset_name=sg.dataset.name)

        expected_outputs = self.expected_outputs(sg)

        # instead of handling, we should probably just exclude this and run again
        if sg.id not in query_result:
            raise ValueError(f'No SV VCFs recorded for {sg.id}')

        lr_sv_vcf: str = query_result[sg.id]
        local_vcf = get_batch().read_input(lr_sv_vcf)

        mod_job = get_batch().new_bash_job(f'Convert {lr_sv_vcf} prior to annotation')
        mod_job.storage('10Gi')
        mod_job.image(config_retrieve(['workflow', 'driver_image']))

        # mandatory argument
        ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
        fasta = get_batch().read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

        # the console entrypoint for the sniffles modifier script has only existed since 1.25.13, requires >=1.25.13
        sex = sg.pedigree.sex.value if sg.pedigree.sex else 0
        mod_job.command(
            'modify_sniffles '
            f'--vcf_in {local_vcf} '
            f'--vcf_out {mod_job.output} '
            f'--new_id {sg.id} '
            f'--fa {fasta} '
            f'--sex {sex} '
            f'--sv',
        )

        # block-gzip and index that result
        tabix_job = get_batch().new_job(f'BGZipping and Indexing for {sg.id}', {'tool': 'bcftools'})
        tabix_job.declare_resource_group(vcf_out={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        tabix_job.image(image=image_path('bcftools'))
        tabix_job.storage('10Gi')
        tabix_job.command(f'bcftools view {mod_job.output} | bgzip -c > {tabix_job.vcf_out["vcf.bgz"]}')
        tabix_job.command(f'tabix {tabix_job.vcf_out["vcf.bgz"]}')

        # write from temp storage into GCP
        get_batch().write_output(tabix_job.vcf_out, str(expected_outputs['vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(target=sg, jobs=[mod_job, tabix_job], data=expected_outputs)


@stage(required_stages=ReFormatPacBioSVs)
class MergeLongReadSVs(MultiCohortStage):
    """
    find all the amended VCFs, and do a naive merge into one huge VCF
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'vcf': self.prefix / 'merged_reformatted_svs.vcf.bgz',
            'index': self.prefix / 'merged_reformatted_svs.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        reuse the gVCF method which does the same? FastCombineVcfs. Nah, not yet - that contains other modifications
        """

        outputs = self.expected_outputs(multicohort)

        # do a bcftools merge on all input files
        modified_vcfs = inputs.as_dict_by_target(ReFormatPacBioSVs)

        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in [str(modified_vcfs[sgid]['vcf']) for sgid in multicohort.get_sequencing_group_ids()]
        ]

        merge_job = get_batch().new_job('Merge Long-Read SV calls', attributes={'tool': 'bcftools'})
        merge_job.image(image=image_path('bcftools'))

        # guessing at resource requirements
        merge_job.cpu(4)
        merge_job.memory('16Gi')
        merge_job.storage('50Gi')
        merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # option breakdown:
        # -Oz: bgzip output
        # -o: output file
        # --threads: number of threads to use
        # -m: merge strategy
        # -0: compression level
        merge_job.command(
            f'bcftools merge {" ".join(batch_vcfs)} -Oz -o '
            f'{merge_job.output["vcf.bgz"]} --threads 4 -m none -0',  # type: ignore
        )
        merge_job.command(f'tabix {merge_job.output["vcf.bgz"]}')  # type: ignore

        # write the result out
        get_batch().write_output(merge_job.output, str(outputs['vcf']).removesuffix('.vcf.bgz'))  # type: ignore

        return self.make_outputs(multicohort, data=outputs, jobs=merge_job)


@stage(required_stages=MergeLongReadSVs, analysis_type='sv', analysis_keys=['annotated_vcf'])
class AnnotateLongReadSVs(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'annotated_vcf': self.prefix / 'annotated_long_read_svs.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'annotated_long_read_svs.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        use the communal GATK-SV wrapper to annotate this merged VCF
        """
        expected_out = self.expected_outputs(multicohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        job_or_none = queue_annotate_sv_jobs(
            multicohort=multicohort,
            prefix=self.prefix,
            input_vcf=inputs.as_dict(multicohort, MergeLongReadSVs)['vcf'],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(multicohort, data=expected_out, jobs=job_or_none)


@stage(required_stages=AnnotateLongReadSVs)
class AnnotateLongReadSVsWithStrvctvre(MultiCohortStage):
    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'strvctvre_vcf': self.prefix / 'strvctvre_annotated.vcf.gz',
            'strvctvre_vcf_index': self.prefix / 'strvctvre_annotated.vcf.gz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:

        input_dict = inputs.as_dict(multicohort, AnnotateLongReadSVs)
        input_vcf = get_batch().read_input_group(
            vcf=str(input_dict['annotated_vcf']),
            vcf_index=str(input_dict['annotated_vcf_index']),
        )['vcf']
        expected_out = self.expected_outputs(multicohort)

        strvctvre_job = queue_annotate_strvctvre_job(
            input_vcf,
            str(expected_out['strvctvre_vcf']),
            self.get_job_attrs(),
        )

        return self.make_outputs(multicohort, data=expected_out, jobs=strvctvre_job)


@stage(required_stages=AnnotateLongReadSVsWithStrvctvre)
class AnnotateCohortLRSv(MultiCohortStage):
    """
    First step to transform annotated SV callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Expected to write a matrix table.
        """
        return {'mt': self.prefix / 'cohort_sv.mt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)

        vcf_path = inputs.as_path(target=multicohort, stage=AnnotateLongReadSVsWithStrvctvre, key='strvctvre_vcf')

        job = annotate_cohort_jobs_sv(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            checkpoint_prefix=self.tmp_prefix / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[AnnotateCohortLRSv], analysis_type='sv', analysis_keys=['mt'])
class AnnotateDatasetLRSv(DatasetStage):
    """
    Subset the MT to be this Dataset only
    Then work up all the genotype values
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Expected to generate a matrix table
        """

        return {'mt': (dataset.prefix() / 'mt' / f'LongReadSV-{get_workflow().output_version}-{dataset.name}.mt')}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortLRSv, key='mt')

        outputs = self.expected_outputs(dataset)

        checkpoint_prefix = dataset.tmp_prefix() / dataset.name / 'checkpoints'

        jobs = annotate_dataset_jobs_sv(
            mt_path=mt_path,
            sgids=dataset.get_sequencing_group_ids(),
            out_mt_path=outputs['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage(
    required_stages=[AnnotateDatasetLRSv],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SV'},
)
class MtToEsLrSv(DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-LR-SV-{get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Uses the non-DataProc MT-to-ES conversion script
        """

        # try to generate a password here - we'll find out inside the script anyway, but
        # by that point we'd already have localised the MT, wasting time and money
        try:
            _es_password_string = es_password()
        except PermissionDenied:
            get_logger().warning(f'No permission to access ES password, skipping for {dataset}')
            return self.make_outputs(dataset)
        except KeyError:
            get_logger().warning(f'ES section not in config, skipping for {dataset}')
            return self.make_outputs(dataset)

        outputs = self.expected_outputs(dataset)

        # get the absolute path to the MT
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetLRSv, key='mt'))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_job(f'Generate {index_name} from {mt_path}')
        # set all job attributes in one bash
        job.cpu(4).memory('lowmem').storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))

        # localise the MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

        # run the export from the localised MT - this job writes no new data, just transforms and exports over network
        job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
