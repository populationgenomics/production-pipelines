"""
Stage that performs joint genotyping of GVCFs using GATK.
"""
import logging

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.filetypes import GvcfPath
from cpg_workflows.jobs import joint_genotyping, seqr_loader
from cpg_workflows.workflow import (
    Cohort,
    stage,
    StageInput,
    StageOutput,
    CohortStage,
    WorkflowError,
)
from .genotype import Genotype
from .. import get_batch
from ..resources import joint_calling_scatter_count


@stage(required_stages=Genotype)
class JointGenotyping(CohortStage):
    """
    Joint-calling of GVCFs together.

    Optionally uses traditional joint calling, or the faster vcf combiner approach.
    """

    def expected_outputs(self, cohort: Cohort) -> dict:
        """
        Generate a site-only VCF and a hail mt.
        """
        outputs = {
            # writing into perm location for late debugging
            # convert to str to avoid checking existence
            'tmp_prefix': str(self.prefix / 'tmp'),
            'mt': self.prefix / 'full.mt',
            'siteonly': self.prefix / 'siteonly.vcf.bgz',
        }

        if get_config()['workflow'].get('use_vcf_combiner', False):
            outputs['vds'] = str(self.prefix / 'full.vds')
            outputs['vcf'] = outputs['mt']
        else:
            outputs['vcf'] = to_path(self.prefix / 'full.vcf.gz')
            outputs['siteonly_part_pattern'] = str(
                self.prefix / 'siteonly_parts' / 'part{idx}.vcf.gz'
            ),

        return outputs

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Submit jobs.
        """
        gvcf_by_sgid = {
            sequencing_group.id: GvcfPath(
                inputs.as_path(target=sequencing_group, stage=Genotype, key='gvcf')
            )
            for sequencing_group in cohort.get_sequencing_groups()
        }

        not_found_gvcfs: list[str] = []
        for sgid, gvcf_path in gvcf_by_sgid.items():
            if gvcf_path is None:
                logging.error(f'Joint genotyping: could not find GVCF for {sgid}')
                not_found_gvcfs.append(sgid)
        if not_found_gvcfs:
            raise WorkflowError(
                f'Joint genotyping: could not find {len(not_found_gvcfs)} '
                f'GVCFs, exiting'
            )

        jobs = []
        siteonly_vcf_path = self.expected_outputs(cohort)['siteonly']
        scatter_count = joint_calling_scatter_count(len(cohort.get_sequencing_groups()))

        if get_config()['workflow'].get('use_vcf_combiner', False):
            from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
            from cpg_workflows.large_cohort.combiner import run

            # Combine GVCFs using vcf combiner
            combine_job = dataproc_job(
                    job_name=self.__class__.__name__,
                    function=run,
                    function_path_args=dict(
                        out_vds_path=self.expected_outputs(cohort)['vds'],
                        tmp_prefix=self.tmp_prefix,
                    ),
                    autoscaling_policy=(
                        get_config()['hail']
                        .get('dataproc', {})
                        .get('combiner_autoscaling_policy')
                    ),
                    num_workers=scatter_count,
                    depends_on=inputs.get_jobs(cohort),
                )

            jobs.append(combine_job)

            # Convert VDS to MT and sites-only VCF
            to_mt_job = seqr_loader.vds_to_mt_job(
                b=get_batch(),
                vds_path=self.expected_outputs(cohort)['vds'],
                out_mt_path=self.expected_outputs(cohort)['mt'],
                out_siteonly_vcf_path=siteonly_vcf_path,
                depends_on=[combine_job],
            )
            jobs.append(to_mt_job)

            from cpg_workflows.large_cohort.site_only_vcf import run as site_only_run

            sites_only_job = dataproc_job(
                job_name=self.__class__.__name__,
                function=site_only_run,
                function_path_args=dict(
                    vds_path=self.expected_outputs(cohort)['vds'],
                    out_vcf_path=siteonly_vcf_path,
                    tmp_prefix=self.tmp_prefix,
                ),
                depends_on=[combine_job],
                # hl.export_vcf() uses non-preemptible workers' disk to merge VCF files.
                # 10 samples take 2.3G, 400 samples take 60G, which roughly matches
                # `huge_disk` (also used in the AS-VQSR VCF-gather job)
                worker_boot_disk_size=200,
                secondary_worker_boot_disk_size=200,
            )
            jobs.append(sites_only_job)

        else:
            # Run joint genotyping using GenotypeGVCFs or GnarlyGenotyper.
            vcf_path = self.expected_outputs(cohort)['vcf']
            out_siteonly_vcf_part_paths = [
                to_path(
                    self.expected_outputs(cohort)['siteonly_part_pattern'].format(
                        idx=idx
                    )
                )
                for idx in range(scatter_count)
            ]

            jc_jobs = joint_genotyping.make_joint_genotyping_jobs(
                b=get_batch(),
                out_vcf_path=vcf_path,
                out_siteonly_vcf_path=siteonly_vcf_path,
                tmp_bucket=to_path(self.expected_outputs(cohort)['tmp_prefix']),
                gvcf_by_sgid=gvcf_by_sgid,
                tool=joint_genotyping.JointGenotyperTool.GnarlyGenotyper
                if get_config()['workflow'].get('use_gnarly', False)
                else joint_genotyping.JointGenotyperTool.GenotypeGVCFs,
                intervals_path=get_config()['workflow'].get('intervals_path'),
                job_attrs=self.get_job_attrs(),
                scatter_count=scatter_count,
                out_siteonly_vcf_part_paths=out_siteonly_vcf_part_paths,
            )
            jobs.extend(jc_jobs)

        for job in jobs:
            assert job

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
