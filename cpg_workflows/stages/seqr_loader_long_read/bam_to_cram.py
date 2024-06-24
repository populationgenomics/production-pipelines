"""
Stage that converts a BAM file to a CRAM file.
Intended for use with long-read BAM files from PacBio.
"""
from functools import cache
from cpg_utils import Path
from cpg_utils.config import config_retrieve, image_path, reference_path, try_get_ar_guid, AR_GUID_NAME
from cpg_utils.hail_batch import get_batch
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import bam_to_cram, seqr_loader_long_read
from cpg_workflows.utils import ExpectedResultT
from cpg_workflows.workflow import (
    Cohort,
    CohortStage,
    SequencingGroup,
    SequencingGroupStage,
    StageInput,
    StageOutput,
    stage,
)
from cpg_workflows.stages.gatk_sv.gatk_sv_common import queue_annotate_sv_jobs

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
}
    """
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
    for sg_id_section in analysis_results['data']['project']['sequencingGroups']:
        for analysis in sg_id_section['analyses']:
            if analysis['output'].endswith('.SVs.phased.vcf.gz'):
                return_dict[sg_id_section['id']] = analysis['output']

    return return_dict


def make_long_read_cram_path(sequencing_group: SequencingGroup) -> CramPath:
    """
    Path to a CRAM file. Not checking its existence here.
    """
    path: Path = sequencing_group.dataset.prefix() / 'pacbio' / 'cram' / f'{sequencing_group.id}.cram'
    return CramPath(
        path=path,
        index_path=path.with_suffix('.cram.crai'),
        reference_assembly=reference_path('broad/ref_fasta'),
    )


@stage(analysis_type=config_retrieve(['workflow', 'bam_to_cram_analysis_type'], 'cram'), analysis_keys=['cram'])
class BamToCram(SequencingGroupStage):
    """
    Convert a BAM to a CRAM file.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        if sequencing_group.cram:
            cram_path = sequencing_group.cram
        else:
            cram_path = make_long_read_cram_path(sequencing_group)

        return {'cram': cram_path.path, 'cram.crai': cram_path.index_path}

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Using the existing `bam_to_cram` function from the `jobs` module.
        """
        b = get_batch()
        input_bam = b.read_input_group(bam=str(sequencing_group.alignment_input_by_seq_type.get('genome')))
        job, output_cram = bam_to_cram.bam_to_cram(
            b=get_batch(),
            input_bam=input_bam,
            extra_label='long_read',
            job_attrs=self.get_job_attrs(sequencing_group),
            requested_nthreads=1,
        )
        b.write_output(output_cram, str(self.expected_outputs(sequencing_group)['cram']).removesuffix('.cram'))

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=[job])


@stage(analysis_keys=['vcf'], analysis_type='custom')
class ReFormatPacBioSVs(SequencingGroupStage):
    """
    take each of the long-read SV VCFs, and re-format the contents
    the huge REF strings are just not required
    a symbolic ALT allele (instead of "N") is required for annotation
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        return {
            'vcf': self.prefix / f'{sequencing_group.id}_reformatted_svs.vcf.bgz',
            'index': self.prefix / f'{sequencing_group.id}_reformatted_svs.vcf.bgz.tbi'
        }

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput:
        """
        a python job to change the VCF contents
        - use a python job to modify all VCF lines, in-line
        - use a follow-up job to block-gzip and index
        """

        # find the vcf for this SG
        query_result = query_for_sv_vcfs(dataset_name=sequencing_group.dataset.name)

        expected_outputs = self.expected_outputs(sequencing_group)

        # instead of handling, we should probably just exclude this and run again
        if sequencing_group.id not in query_result:
            raise ValueError(f'No SV VCFs recorded for {sequencing_group.id}')

        lr_sv_vcf: str = query_result[sequencing_group.id]
        local_vcf = get_batch().read_input(lr_sv_vcf)

        py_job = get_batch().new_python_job(f'Convert {lr_sv_vcf} prior to annotation')
        py_job.storage('10Gi')
        py_job.call(seqr_loader_long_read.modify_sniffles_vcf, local_vcf, py_job.output)

        # block-gzip and index that result
        tabix_job = get_batch().new_job(
            name=f'BGZipping and Indexing for {sequencing_group.id}',
            attributes={'tool': 'bcftools'},
        )
        tabix_job.declare_resource_group(vcf_out={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        tabix_job.image(image=image_path('bcftools'))
        tabix_job.storage('10Gi')
        tabix_job.command(f'bcftools view {py_job.output} | bgzip -c > {tabix_job.vcf_out["vcf.bgz"]}')  # type: ignore
        tabix_job.command(f'tabix {tabix_job.vcf_out["vcf.bgz"]}')  # type: ignore

        # write from temp storage into GCP
        get_batch().write_output(tabix_job.vcf_out, str(expected_outputs['vcf']).removesuffix('vcf.bgz'))

        return self.make_outputs(target=sequencing_group, jobs=[py_job, tabix_job], data=expected_outputs)


@stage(required_stages=ReFormatPacBioSVs)
class MergeLongReadSVs(CohortStage):
    """
    find all the amended VCFs, and do a naive merge into one huge VCF
    """
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        return {
            'vcf': self.prefix / 'merged_reformatted_svs.vcf.bgz',
            'index': self.prefix / 'merged_reformatted_svs.vcf.bgz.tbi'
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        reuse the gVCF method which does the same? FastCombineVcfs. Nah, not yet - that contains other modifications
        """

        outputs = self.expected_outputs(cohort)

        # do a slapdash bcftools merge on all input files...
        modified_vcfs = inputs.as_dict_by_target(ReFormatPacBioSVs)

        batch_vcfs = [
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz']
            for each_vcf in [str(modified_vcfs[sgid]['vcf']) for sgid in cohort.get_sequencing_group_ids()]
        ]

        # just gonna create a job and commands here... this could be in a jobs file, but... shrug
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
            f'{merge_job.output["vcf.bgz"]} --threads 4 -m all -0'  # type: ignore
        )
        merge_job.command(f'tabix {merge_job.output["vcf.bgz"]}')  # type: ignore

        # write the result out
        get_batch().write_output(merge_job.output, outputs['vcf'].removesuffix('vcf.bgz'))

        return self.make_outputs(cohort, data=outputs, jobs=merge_job)


@stage(required_stages=MergeLongReadSVs, analysis_type='custom', analysis_keys=['annotated_vcf'])
class AnnotateLongReadSVs(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'annotated_vcf': self.prefix / 'annotated_long_read_svs.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'annotated_long_read_svs.vcf.bgz.tbi'
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        use the communal GATK-SV wrapper to annotate this merged VCF
        """
        expected_out = self.expected_outputs(cohort)

        billing_labels = {'stage': self.name.lower(), AR_GUID_NAME: try_get_ar_guid()}

        job_or_none = queue_annotate_sv_jobs(
            cohort=cohort,
            cohort_prefix=self.prefix,
            input_vcf=inputs.as_dict(cohort, MergeLongReadSVs)['vcf'],
            outputs=expected_out,
            labels=billing_labels,
        )
        return self.make_outputs(cohort, data=expected_out, jobs=job_or_none)


