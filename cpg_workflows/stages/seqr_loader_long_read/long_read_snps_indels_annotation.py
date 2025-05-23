"""
Workflow for finding long-read SNPs_Indels files, updating file contents, merging, and annotating the results
"""

from functools import cache

from google.api_core.exceptions import PermissionDenied

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch
from cpg_workflows.jobs.seqr_loader import annotate_dataset_jobs
from cpg_workflows.jobs.seqr_loader_snps_indels import (
    annotate_cohort_jobs_snps_indels,
    split_merged_vcf_and_get_sitesonly_vcfs_for_vep,
)
from cpg_workflows.jobs.vep import add_vep_jobs
from cpg_workflows.resources import joint_calling_scatter_count
from cpg_workflows.stages.seqr_loader import es_password
from cpg_workflows.targets import Dataset, MultiCohort, SequencingGroup
from cpg_workflows.utils import get_logger, tshirt_mt_sizing
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
        meta
        output
        outputs
      }
    }
  }
}""",
)

LRS_IDS_QUERY = gql(
    """
    query MyQuery($dataset: String!) {
    project(name: $dataset) {
      sequencingGroups(technology: {eq: "long-read"}, type: {eq: "genome"}) {
          id
          sample {
            externalId
            meta
            participant {
              externalId
              reportedSex
            }
          }
        }
      }
    }
    """,
)


@cache
def query_for_snps_indels_vcfs(
    dataset_name: str,
    pipeface_versions: tuple[str] | None,
    snps_indels_callers: tuple[str] | None,
) -> dict[str, dict]:
    """
    query metamist for the PacBio SNPs_Indels VCFs
    return a dictionary of each CPG ID and its corresponding VCF
    this is cached - used in a SequencingGroupStage, but we only want to query for it once instead of once/SG

    Args:
        dataset_name (str): the name of the dataset
        pipeface_versions (tuple[str] | None): a tuple of allowed pipeface versions (hashable, so can be cached)

    Returns:
        a dictionary of the SG IDs and their phased SNPs Indels VCF
    """
    single_sample_vcfs: dict[str, dict] = {}
    joint_called_vcfs: dict[str, dict] = {}
    analysis_results = query(VCF_QUERY, variables={'dataset': dataset_name})
    for sg in analysis_results['project']['sequencingGroups']:
        for analysis in sg['analyses']:
            if pipeface_versions and analysis['meta'].get('pipeface_version', 'unkown') not in pipeface_versions:
                get_logger().info(
                    f'Skipping {analysis["output"]} for {sg["id"]} as it is not an allowed pipeface version: {pipeface_versions}',
                )
                continue
            if snps_indels_callers and analysis['meta'].get('caller', 'unkown') not in snps_indels_callers:
                get_logger().info(
                    f'Skipping {analysis["output"]} for {sg["id"]} as it is not an allowed caller: {snps_indels_callers}',
                )
                continue
            if not analysis['output'].endswith('snp_indel.phased.vcf.gz'):
                continue
            if analysis['meta'].get('joint_called', False):
                joint_called_vcfs[sg['id']] = {
                    'output': analysis['output'],
                    'meta': analysis['meta'],
                }
            else:
                single_sample_vcfs[sg['id']] = {
                    'output': analysis['output'],
                    'meta': analysis['meta'],
                }

    # Prefer the joint-called VCFs over the single-sample VCFs
    sg_vcfs = {}
    for sg_id, single_sample_vcf in single_sample_vcfs.items():
        if sg_id not in joint_called_vcfs:
            sg_vcfs[sg_id] = single_sample_vcf
            continue
        sg_vcfs[sg_id] = joint_called_vcfs[sg_id]
    for sg_id, joint_called_vcf in joint_called_vcfs.items():
        if sg_id not in sg_vcfs:
            sg_vcfs[sg_id] = joint_called_vcf

    # Remove the parents entries if their family has a joint-called VCF
    sgs_to_skip = find_sgs_to_skip(sg_vcfs)
    return_dict = {}
    for sg_id, vcf_analysis in sg_vcfs.items():
        if sg_id in sgs_to_skip:
            get_logger().info(f'Skipping {sg_id} as it is a parent in a joint-called VCF')
            continue
        return_dict[sg_id] = vcf_analysis

    return return_dict


def find_sgs_to_skip(sg_vcf_dict: dict[str, dict]) -> set[str]:
    """
    Find the SGs to skip in the reformatting stage
    These are the parents if the family has been joint-called
    """
    sgs_to_skip = set()
    joint_called_families = set()
    for sg_id, vcf_analysis in sg_vcf_dict.items():
        analysis_meta = vcf_analysis['meta']
        if analysis_meta.get('joint_called', False):
            joint_called_families.add(analysis_meta.get('family_id', ''))
    for sg_id, vcf_analysis in sg_vcf_dict.items():
        analysis_meta = vcf_analysis['meta']
        if analysis_meta.get('family_id', '') in joint_called_families and not analysis_meta.get('joint_called', False):
            sgs_to_skip.add(sg_id)
    return sgs_to_skip


def query_for_lrs_sg_id_and_sex_mapping(datasets: list[str]):
    """
    Query metamist for the LRS ID corresponding to each sequencing group ID, and to its participant's sex
    """
    lrs_sgid_mapping = {}
    lrs_id_sex_mapping = {}
    for dataset in datasets:
        if config_retrieve(['workflow', 'access_level']) == 'test' and not dataset.endswith('-test'):
            dataset = dataset + '-test'
        query_results = query(LRS_IDS_QUERY, variables={'dataset': dataset})
        for sg in query_results['project']['sequencingGroups']:
            sample = sg['sample']
            participant = sample['participant']
            lrs_id = sample['meta'].get('lrs_id', None)
            if not lrs_id:
                get_logger().warning(
                    f'{dataset} :: No LRS ID found for {participant["externalId"]} - {sample["externalId"]}',
                )
                continue
            lrs_sgid_mapping[lrs_id] = sg['id']
            lrs_id_sex_mapping[lrs_id] = participant['reportedSex']
    return lrs_sgid_mapping, lrs_id_sex_mapping


@stage
class WriteLRSIDtoSGIDMappingFile(MultiCohortStage):
    """
    Write the LRS ID to SG ID mapping to a file
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'lrs_sgid_mapping': self.prefix / 'lrs_sgid_mapping.txt',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Write the LRS ID to SG ID mapping to a file
        This is used by bcftools reheader to update the sample IDs in the VCFs
        """
        output = self.expected_outputs(multicohort)

        lrs_sgid_mapping, _ = query_for_lrs_sg_id_and_sex_mapping([d.name for d in multicohort.get_datasets()])
        mapping_file_path = self.prefix / 'lrs_sgid_mapping.txt'
        get_logger().info(f'Writing LRS ID to SG ID mapping to {mapping_file_path}')
        with mapping_file_path.open('w') as f:
            for lrs_id, sg_id in lrs_sgid_mapping.items():
                f.write(f'{lrs_id} {sg_id}\n')

        return self.make_outputs(multicohort, data=output)


@stage(required_stages=[WriteLRSIDtoSGIDMappingFile])
class ReFormatPacBioSNPsIndels(SequencingGroupStage):
    """
    take each of the long-read SNPs Indels VCFs, and re-format the contents
    the huge REF strings are just not required
    a symbolic ALT allele (instead of "N") is required for annotation
    adds a unique ID to each record, so that GATK-SV can do its weird sorting (is this necessary for SNPs?)
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        sgid_prefix = sequencing_group.dataset.prefix() / 'pacbio' / 'modified_vcfs'
        return {
            'vcf': sgid_prefix / f'{sequencing_group.id}_reformatted_renamed_lr_snps_indels.vcf.bgz',
            'index': sgid_prefix / f'{sequencing_group.id}_reformatted_renamed_lr_snps_indels.vcf.bgz.tbi',
        }

    def queue_jobs(self, sg: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        a python job to change the VCF contents
        - use a python job to modify all VCF lines, in-line
        - a bcftools job to reheader the VCF with the replacement sample IDs, and then sort
        - use a follow-up job to block-gzip
        - this is skipped for the parents in trio joint-called VCFs
        """
        # find the VCFs for this dataset
        dataset_name = (
            sg.dataset.name + '-test' if config_retrieve(['workflow', 'access_level']) == 'test' else sg.dataset.name
        )
        pipeface_versions = config_retrieve(
            ['workflow', 'long_read_vcf_annotation', 'pipeface_versions'],
            default=None,
        )
        if pipeface_versions:
            pipeface_versions = tuple(pipeface_versions)
        snps_indels_callers = config_retrieve(
            ['workflow', 'long_read_vcf_annotation', 'snps_indels_callers'],
            default=None,
        )
        if snps_indels_callers:
            snps_indels_callers = tuple(snps_indels_callers)
        sg_vcfs = query_for_snps_indels_vcfs(
            dataset_name=dataset_name,
            pipeface_versions=pipeface_versions,
            snps_indels_callers=snps_indels_callers,
        )
        if sg.id not in sg_vcfs:
            return None

        expected_outputs = self.expected_outputs(sg)

        lr_snps_indels_vcf: str = sg_vcfs[sg.id]['output']
        local_vcf = get_batch().read_input(lr_snps_indels_vcf)

        lrs_sg_id_map = inputs.as_path(get_multicohort(), WriteLRSIDtoSGIDMappingFile, 'lrs_sgid_mapping')
        local_id_mapping = get_batch().read_input(lrs_sg_id_map)

        joint_called = sg_vcfs[sg.id]['meta'].get('joint_called', False)

        # mandatory argument
        ref_fasta = config_retrieve(['workflow', 'ref_fasta'])
        fasta = get_batch().read_input_group(**{'fa': ref_fasta, 'fa.fai': f'{ref_fasta}.fai'})['fa']

        # normalise, reheader, then block-gzip and index that result
        tabix_job = get_batch().new_job(
            f'Normalising, Reheadering, BGZipping, and Indexing for {sg.id}: {"joint-called " if joint_called else ""}{lr_snps_indels_vcf}',
            {'tool': 'bcftools'},
        )
        tabix_job.declare_resource_group(vcf_out={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        tabix_job.image(image=image_path('bcftools'))
        tabix_job.storage('10Gi')
        tabix_job.command(
            f'bcftools view -Ov {local_vcf} | bcftools reheader --samples {local_id_mapping} -o {tabix_job.reheadered}',
        )
        tabix_job.command(
            f'bcftools norm -m -any -f {fasta} -c s {tabix_job.reheadered} | bcftools sort | bgzip -c > {tabix_job.vcf_out["vcf.bgz"]}',
        )
        tabix_job.command(f'tabix {tabix_job.vcf_out["vcf.bgz"]}')

        # write from temp storage into GCP
        get_batch().write_output(tabix_job.vcf_out, str(expected_outputs['vcf']).removesuffix('.vcf.bgz'))

        return self.make_outputs(target=sg, jobs=[tabix_job], data=expected_outputs)


@stage(required_stages=ReFormatPacBioSNPsIndels)
class MergeLongReadSNPsIndels(MultiCohortStage):
    """
    Find all the amended VCFs, and do a naive merge into one huge VCF
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        return {
            'vcf': self.prefix / 'merged_reformatted_snps_indels.vcf.bgz',
            'index': self.prefix / 'merged_reformatted_snps_indels.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Use bcftools to merge all the VCFs, and then fill in the tags (requires bcftools 1.12+)
        """

        outputs = self.expected_outputs(multicohort)

        # do a bcftools merge on all input files
        modified_vcfs = inputs.as_dict_by_target(ReFormatPacBioSNPsIndels)
        if len(modified_vcfs) == 1:
            get_logger().info('Only one VCF found, skipping merge')
            return None

        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in [
                str(modified_vcfs.get(sgid, {}).get('vcf'))
                for sgid in multicohort.get_sequencing_group_ids()
                if sgid in modified_vcfs
            ]
        ]

        merge_job = get_batch().new_job('Merge Long-Read SNPs Indels calls', attributes={'tool': 'bcftools'})
        merge_job.image(image=image_path('bcftools_120'))

        # guessing at resource requirements
        merge_job.cpu(config_retrieve(['resource_overrides', 'merge_vcfs', 'ncpu'], 4))
        merge_job.memory(str(config_retrieve(['resource_overrides', 'merge_vcfs', 'mem_gb'], '16')) + 'Gi')
        merge_job.storage(str(config_retrieve(['resource_overrides', 'merge_vcfs', 'storage_gb'], '50')) + 'Gi')
        merge_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # option breakdown:
        # -Oz: bgzip output
        # -o: output file
        # --threads: number of threads to use
        # -m: merge strategy
        # -0: compression level
        merge_job.command(
            f'bcftools merge {" ".join(batch_vcfs)} -Oz -o temp.vcf.bgz --threads 4 -m none -0',  # type: ignore
        )
        merge_job.command(
            f'bcftools +fill-tags temp.vcf.bgz -Oz -o {merge_job.output["vcf.bgz"]} --write-index=tbi -- -t AF,AN,AC',
        )

        # write the result out
        get_batch().write_output(merge_job.output, str(outputs['vcf']).removesuffix('.vcf.bgz'))  # type: ignore

        return self.make_outputs(multicohort, data=outputs, jobs=merge_job)


@stage(required_stages=[ReFormatPacBioSNPsIndels, MergeLongReadSNPsIndels])
class SiteOnlyVCFs(MultiCohortStage):
    """
    Get the site-only VCFs from the merged VCF
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict:
        """
        Generate site-only VCFs from the merged VCF.
        """
        return {
            'siteonly': to_path(self.prefix / 'siteonly.vcf.gz'),
            'siteonly_part_pattern': str(self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'),
        }

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        jobs = []
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))
        out_siteonly_vcf_part_paths = [
            to_path(outputs['siteonly_part_pattern'].format(idx=idx)) for idx in range(scatter_count)
        ]

        intervals_path = None
        if config_retrieve(['workflow', 'intervals_path'], default=None):
            intervals_path = to_path(config_retrieve(['workflow', 'intervals_path']))

        exclude_intervals_path = None
        if config_retrieve(['workflow', 'exclude_intervals_path'], default=None):
            exclude_intervals_path = to_path(config_retrieve(['workflow', 'exclude_intervals_path']))

        if len(inputs.as_dict_by_target(ReFormatPacBioSNPsIndels)) == 1:
            merged_vcf_path = inputs.as_path(multicohort, ReFormatPacBioSNPsIndels, 'vcf')
        else:
            merged_vcf_path = inputs.as_path(multicohort, MergeLongReadSNPsIndels, 'vcf')

        vcf_jobs = split_merged_vcf_and_get_sitesonly_vcfs_for_vep(
            b=get_batch(),
            scatter_count=scatter_count,
            merged_vcf_path=merged_vcf_path,
            tmp_bucket=self.tmp_prefix / 'tmp',
            out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
            intervals_path=intervals_path,
            exclude_intervals_path=exclude_intervals_path,
            job_attrs=self.get_job_attrs(),
        )
        jobs.extend(vcf_jobs)

        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage(required_stages=[SiteOnlyVCFs])
class VepLRS(MultiCohortStage):
    """
    Run VEP on the long-read site-only VCFs and write out a Hail table.
    """

    def expected_outputs(self, multicohort: MultiCohort):
        """
        Expected to write a hail table.
        """
        return {'ht': self.prefix / 'long_read' / 'vep.ht'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        outputs = self.expected_outputs(multicohort)
        scatter_count = joint_calling_scatter_count(len(multicohort.get_sequencing_groups()))
        input_siteonly_vcf_part_paths = [
            to_path(
                inputs.as_str(
                    stage=SiteOnlyVCFs,
                    target=multicohort,
                    key='siteonly_part_pattern',
                ).format(idx=idx),
            )
            for idx in range(scatter_count)
        ]

        jobs = add_vep_jobs(
            get_batch(),
            input_vcfs=input_siteonly_vcf_part_paths,
            out_path=outputs['ht'],
            tmp_prefix=self.tmp_prefix / 'tmp',
            job_attrs=self.get_job_attrs(),
            scatter_count=scatter_count,
        )
        return self.make_outputs(multicohort, outputs, jobs)


@stage(required_stages=[ReFormatPacBioSNPsIndels, MergeLongReadSNPsIndels, VepLRS])
class AnnotateCohortLRSNPsIndels(MultiCohortStage):
    """
    First step to transform annotated SNPs Indels callset data into a seqr ready format
    """

    def expected_outputs(self, multicohort: MultiCohort) -> dict[str, Path]:
        """
        Expected to write a matrix table.
        """
        return {'mt': self.prefix / 'cohort_snps_indels.mt'}

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput | None:
        """
        queue job(s) to rearrange the annotations prior to Seqr transformation
        """
        outputs = self.expected_outputs(multicohort)

        if len(inputs.as_dict_by_target(ReFormatPacBioSNPsIndels)) == 1:
            vcf_path = inputs.as_path(multicohort, ReFormatPacBioSNPsIndels, 'vcf')
        else:
            vcf_path = inputs.as_path(target=multicohort, stage=MergeLongReadSNPsIndels, key='vcf')

        vep_ht_path = inputs.as_path(target=multicohort, stage=VepLRS, key='ht')

        job = annotate_cohort_jobs_snps_indels(
            vcf_path=vcf_path,
            out_mt_path=outputs['mt'],
            vep_ht_path=vep_ht_path,
            checkpoint_prefix=self.tmp_prefix / 'checkpoints',
            job_attrs=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=outputs, jobs=job)


@stage(required_stages=[AnnotateCohortLRSNPsIndels], analysis_type='custom', analysis_keys=['mt'])
class AnnotateDatasetLRSNPsIndels(DatasetStage):
    """
    Subset the MT to be this Dataset only
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """
        Expected to generate a matrix table
        """

        return {
            'mt': (dataset.prefix() / 'mt' / f'LongReadSNPsIndels-{get_workflow().output_version}-{dataset.name}.mt'),
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Subsets the whole MT to this cohort only
        Then brings a range of genotype data into row annotations

        Args:
            dataset (Dataset): SGIDs specific to this dataset/project
            inputs ():
        """
        mt_path = inputs.as_path(target=get_multicohort(), stage=AnnotateCohortLRSNPsIndels, key='mt')

        outputs = self.expected_outputs(dataset)

        checkpoint_prefix = dataset.tmp_prefix() / dataset.name / 'checkpoints'

        jobs = annotate_dataset_jobs(
            mt_path=mt_path,
            sequencing_group_ids=dataset.get_sequencing_group_ids(),
            out_mt_path=outputs['mt'],
            tmp_prefix=checkpoint_prefix,
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=outputs, jobs=jobs)


@stage(
    required_stages=[AnnotateDatasetLRSNPsIndels],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'SNV_INDEL'},
)
class MtToEsLrSNPsIndels(DatasetStage):
    """
    Create a Seqr index
    https://github.com/populationgenomics/metamist/issues/539
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        index_name = f'{dataset.name}-{sequencing_type}-LR-SNPsIndels-{get_workflow().run_timestamp}'.lower()
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
        mt_path = str(inputs.as_path(target=dataset, stage=AnnotateDatasetLRSNPsIndels, key='mt'))
        # and just the name, used after localisation
        mt_name = mt_path.split('/')[-1]

        # get the expected outputs as Strings
        index_name = str(outputs['index_name'])
        flag_name = str(outputs['done_flag'])

        job = get_batch().new_job(f'Generate {index_name} from {mt_path}')
        req_storage = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )
        # set all job attributes in one bash
        job.cpu(4).memory('lowmem').storage(f'{req_storage}Gi').image(config_retrieve(['workflow', 'driver_image']))

        # localise the MT
        job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

        # run the export from the localised MT - this job writes no new data, just transforms and exports over network
        job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {index_name} --flag {flag_name}')

        return self.make_outputs(dataset, data=outputs['index_name'], jobs=job)
