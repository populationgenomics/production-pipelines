#!/usr/bin/env python3

"""
Driver for loading data into SEQR for the CPG. See the README for more information.

- 2021/04/16 Michael Franklin and Vlad Savelyev
"""

import logging
import time
from enum import Enum
from os.path import join, dirname, abspath, splitext, basename
from typing import Optional, List, Tuple, Set, Dict, Collection
import pandas as pd
import click
import hailtop.batch as hb

from cpg_production_pipelines.vqsr import make_vqsr_jobs
from cpg_production_pipelines import utils, resources
from cpg_production_pipelines.jobs import align, split_intervals, haplotype_caller
from cpg_production_pipelines.pipeline import Namespace, Pipeline
from cpg_production_pipelines.smdb import SMDB, parse_reads_from_metadata

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Stage(Enum):
    INPUT = 1
    CRAM = 2
    GVCF = 3
    JOINT_CALLING = 4
    ANNOTATE = 5
    LOAD_TO_ES = 6


@click.command()
@click.option(
    '-n',
    '--namespace',
    'output_namespace',
    type=click.Choice([n.lower() for n in Namespace.__members__]),
    callback=lambda c, p, v: getattr(Namespace, v) if v else None,
    help='The bucket namespace to write the results to',
)
@click.option(
    '--analysis-project',
    'analysis_project',
    default='seqr',
    help='SM project name to write the intermediate/joint-calling analysis entries to',
)
@click.option(
    '--input-project',
    'input_projects',
    multiple=True,
    required=True,
    help='Only read samples that belong to the project(s). Can be set multiple times.',
)
@click.option(
    '--output-project',
    'output_projects',
    multiple=True,
    help='Only create ES indicies for the project(s). Can be set multiple times. '
    'Defaults to --input-projects. The name of the ES index will be suffixed '
    'with the dataset version (set by --version)',
)
@click.option(
    '--start-from-stage',
    'start_from_stage',
    type=click.Choice([n.lower() for n in Stage.__members__]),
    callback=lambda c, p, v: getattr(Stage, v) if v else None,
    help='Only pick results from the previous stages if they exist. '
    'If not, skip such samples',
)
@click.option(
    '--end-with-stage',
    'end_with_stage',
    type=click.Choice([n.lower() for n in Stage.__members__]),
    callback=lambda c, p, v: getattr(Stage, v) if v else None,
    help='Finish the pipeline after this stage',
)
@click.option(
    '--skip-sample',
    '-S',
    'skip_samples',
    multiple=True,
    help='Don\'t process specified samples. Can be set multiple times.',
)
@click.option(
    '--output-version',
    'output_version',
    type=str,
    default='v0',
    help='Suffix the outputs with this version tag. Useful for testing',
)
@click.option('--keep-scratch', 'keep_scratch', is_flag=True)
@click.option(
    '--overwrite/--reuse',
    'overwrite',
    is_flag=True,
    help='if an intermediate or a final file exists, skip running the code '
    'that generates it.',
)
@click.option('--dry-run', 'dry_run', is_flag=True)
@click.option(
    '--make-checkpoints',
    'make_checkpoints',
    is_flag=True,
    help='Create checkpoints for intermediate Hail data',
)
@click.option(
    '--skip-ped-checks',
    'skip_ped_checks',
    is_flag=True,
    help='Skip checking provided sex and pedigree against the inferred one',
)
@click.option('--vep-block-size', 'vep_block_size', type=click.INT)
@click.option(
    '--hc-shards-num',
    'hc_shards_num',
    type=click.INT,
    default=utils.NUMBER_OF_HAPLOTYPE_CALLER_INTERVALS,
    help='Number of intervals to devide the genome for gatk HaplotypeCaller',
)
@click.option(
    '--use-gnarly/--no-use-gnarly',
    'use_gnarly',
    default=False,
    is_flag=True,
    help='Use GnarlyGenotyper instead of GenotypeGVCFs',
)
@click.option(
    '--use-as-vqsr/--no-use-as-vqsr',
    'use_as_vqsr',
    default=True,
    is_flag=True,
    help='Use allele-specific annotations for VQSR',
)
@click.option(
    '--check-inputs-existence/--skip-check-inputs-existence',
    'check_inputs_existence',
    default=True,
    is_flag=True,
)
@click.option(
    '--update-smdb/--skip-update-smdb',
    'update_smdb',
    default=True,
    is_flag=True,
)
def main(
    output_namespace: Namespace,
    analysis_project: str,
    input_projects: Collection[str],
    output_projects: Optional[Collection[str]],
    start_from_stage: Optional[Stage],
    end_with_stage: Optional[Stage],
    skip_samples: Collection[str],
    output_version: str,
    keep_scratch: bool,
    overwrite: bool,
    dry_run: bool,
    skip_ped_checks: bool,  # pylint: disable=unused-argument
    vep_block_size: Optional[int],  # pylint: disable=unused-argument
    hc_shards_num: int,
    use_gnarly: bool,
    use_as_vqsr: bool,
    check_inputs_existence: bool,
    update_smdb: bool,
):  # pylint: disable=missing-function-docstring
    # Determine bucket paths

    assert input_projects
    if output_projects:
        if not all(op in input_projects for op in output_projects):
            logger.critical(
                'All output projects must be contained within '
                'the specified input projects'
            )

    if output_namespace != Namespace.MAIN:
        analysis_project = f'{analysis_project}-test'
        input_projects = [f'{p}-test' for p in input_projects]
        output_projects = [f'{p}-test' for p in output_projects]

    pipeline = SeqrLoaderPipeline(
        analysis_project=analysis_project,
        name='seqr_loader',
        output_version=output_version,
        namespace=output_namespace,
        keep_scratch=keep_scratch,
        title=(
            f'Seqr loading. '
            f'{", ".join(input_projects)} -> '
            f'{", ".join(output_projects)}, '
            f'v{output_version}'
        ),
        do_update_analyses=update_smdb,
        do_check_existence=check_inputs_existence,
    )
    
    pipeline.run(dry_run)


class SeqrLoaderPipeline(Pipeline):
    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.fingerprints_bucket = f'{self.analysis_bucket}/fingerprints'

        # pipeline=pipeline,
        # smdb=smdb,
        # overwrite=overwrite,
        # input_projects=input_projects,
        # output_projects=output_projects or input_projects,
        # vep_block_size=vep_block_size,
        # start_from_stage=start_from_stage,
        # end_with_stage=end_with_stage,
        # skip_samples=skip_samples,
        # use_gnarly=use_gnarly,
        # use_as_vqsr=use_as_vqsr,
        # hc_shards_num=hc_shards_num,
        # check_inputs_existence=check_inputs_existence,

    def add_jobs(
        self,
        input_projects: List[str],
        skip_samples: List[str],
        first_stage: Stage = list(Stage.__members__.keys())[0],
        last_stage: Stage = list(Stage.__members__.keys())[-1],
    ):
        self._prepare_inputs(input_projects, skip_samples)

        if last_stage == Stage.INPUT:
            logger.info(
                f'Latest stage is {last_stage.name}, stopping the pipeline here.'
            )
            return

        if first_stage <= Stage.CRAM:
            self._cram()

    def _cram(
        self,
    ):
        cram_jobs = []
        for proj, samples in self.samples_by_project.items():
            logger.info(f'Submitting CRAMs for project {proj}')
            proj_bucket = f'gs://cpg-{proj}-{self.output_suf}'
            sample_ids = [s['id'] for s in samples]

            cram_analysis_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='cram',
            )
            seq_info_by_sid = self.db.find_seq_info_by_sid(sample_ids)

            for s in samples:
                logger.info(f'Project {proj}. Processing CRAM for {s["id"]}')
                expected_cram_path = f'{proj_bucket}/cram/{s["id"]}.cram'
                found_cram_path = self.db.process_existing_analysis(
                    sample_ids=[s['id']],
                    completed_analysis=cram_analysis_per_sid.get(s['id']),
                    analysis_type='cram',
                    analysis_sample_ids=[s['id']],
                    expected_output_fpath=expected_cram_path,
                )
                seq_info = seq_info_by_sid[s['id']]
                alignment_input = self.db.parse_reads_from_metadata(seq_info['meta'])
                if not alignment_input:
                    logger.critical(f'Could not find read data for sample {s["id"]}')
                    continue
                cram_job = align.bwa(
                    b=self.b,
                    alignment_input=alignment_input,
                    output_path=expected_cram_path,
                    sample_name=s['id'],
                    project_name=proj,
                )
                cram_jobs.append(cram_job)
                found_cram_path = expected_cram_path
        return cram_jobs

    def _gvcf(self, hc_shards_num):
        # after dropping samples with incorrect metadata, missing inputs, etc
        good_samples: List[Dict] = []
        hc_intervals = None
        gvcf_jobs = []
        gvcf_by_sid: Dict[str, str] = dict()
        for proj, samples in self.samples_by_project.items():
            logger.info(f'Processing project {proj}')
            proj_bucket = f'gs://cpg-{proj}-{self.output_suf}'
            sample_ids = [s['id'] for s in samples]

            gvcf_analysis_per_sid = self.db.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type='gvcf',
            )

            for s in samples:
                logger.info(f'Project {proj}. Processing GVCF {s["id"]}')
                expected_gvcf_path = f'{proj_bucket}/gvcf/{s["id"]}.g.vcf.gz'
                found_gvcf_path = self.db.process_existing_analysis(
                    sample_ids=[s['id']],
                    completed_analysis=gvcf_analysis_per_sid.get(s['id']),
                    analysis_type='gvcf',
                    analysis_sample_ids=[s['id']],
                    expected_output_fpath=expected_gvcf_path,
                )
                if hc_intervals is None and hc_shards_num > 1:
                    hc_intervals = split_intervals.intervals(
                        b=self.b,
                        scatter_count=hc_shards_num,
                        ref_fasta=resources.REF_FASTA,
                    )
                gvcf_j = haplotype_caller.produce_gvcf(
                    b=self.b,
                    output_path=expected_gvcf_path,
                    sample_name=s['id'],
                    project_name=proj,
                    cram_path=found_cram_path,
                    crai_path=found_cram_path + '.crai',
                    intervals=hc_intervals,
                    number_of_intervals=hc_shards_num,
                    tmp_bucket=self.tmp_bucket,
                    overwrite=overwrite,
                    depends_on=[cram_job] if cram_job else [],
                    smdb=smdb,
                )
                gvcf_jobs.append(gvcf_j)
                found_gvcf_path = expected_gvcf_path
                gvcf_by_sid[s['id']] = found_gvcf_path
                good_samples.append(s)
        return good_samples

    def _joint_calling(self):
        pass

    def _vqsr(self):
        pass

    def _annotation(self):
        pass

    def _create_index(self):
        pass


def _add_jobs(  # pylint: disable=too-many-statements
    pipeline: Pipeline,
    smdb: SMDB,
    overwrite: bool,
    input_projects: Collection[str],
    output_projects: Collection[str],
    vep_block_size: Optional[int],
    start_from_stage: Optional[str],
    end_with_stage: Optional[str],
    skip_samples: Collection[str],
    use_gnarly: bool,
    use_as_vqsr: bool,
    hc_shards_num: int,
    check_inputs_existence: bool,
) -> Optional[hb.Batch]:

    if end_with_stage == 'gvcf':
        logger.info(f'Latest stage is {end_with_stage}, stopping the pipeline here.')
        return

    if not good_samples:
        logger.info('No samples left to joint-call')
        return

    # Is there a complete joint-calling analysis for the requested set of samples?
    sample_ids = list(set(s['id'] for s in good_samples))
    samples_hash = utils.hash_sample_ids(sample_ids)
    expected_jc_vcf_path = f'{pipeline.tmp_bucket}/joint_calling/{samples_hash}.vcf.gz'
    skip_jc_stage = start_from_stage is not None and start_from_stage not in [
        'cram',
        'gvcf',
        'joint_calling',
    ]
    found_jc_vcf_path = smdb.process_existing_analysis(
        sample_ids=sample_ids,
        completed_analysis=smdb.find_joint_calling_analysis(sample_ids),
        analysis_type='joint-calling',
        analysis_sample_ids=sample_ids,
        expected_output_fpath=expected_jc_vcf_path,
        skip_stage=skip_jc_stage,
    )
    if skip_jc_stage:
        if not found_jc_vcf_path:
            return None
        jc_job = None
    else:
        jc_job = _make_joint_genotype_jobs(
            b=pipeline.b,
            output_path=expected_jc_vcf_path,
            samples=good_samples,
            genomicsdb_bucket=f'{pipeline.analysis_bucket}/genomicsdbs',
            tmp_bucket=pipeline.tmp_bucket,
            gvcf_by_sid=gvcf_by_sid,
            local_tmp_dir=pipeline.local_tmp_dir,
            overwrite=overwrite,
            depends_on=gvcf_jobs,
            analysis_project=pipeline.analysis_project,
            use_gnarly=use_gnarly,
            use_as_vqsr=use_as_vqsr,
        )
        found_jc_vcf_path = expected_jc_vcf_path

    if end_with_stage == 'joint_calling':
        logger.info(f'Latest stage is {end_with_stage}, stopping the pipeline here.')
        return b

    for project in output_projects:
        annotated_mt_path = f'{analysis_bucket}/mt/{project}.mt'

        skip_anno_stage = start_from_stage is not None and start_from_stage not in [
            'cram',
            'gvcf',
            'joint_calling',
            'annotate',
        ]
        if skip_anno_stage:
            annotate_job = None
        else:
            anno_tmp_bucket = f'{tmp_bucket}/annotation/{project}'
            if utils.can_reuse(annotated_mt_path, overwrite):
                annotate_job = b.new_job(f'{project}: annotate [reuse]')
            else:
                sample_map_bucket_path = f'{anno_tmp_bucket}/external_id_map.tsv'
                sample_map_local_fpath = join(
                    local_tmp_dir, basename(sample_map_bucket_path)
                )
                with open(sample_map_local_fpath, 'w') as f:
                    f.write('\t'.join(['s', 'seqr_id']) + '\n')
                    for s in good_samples:
                        f.write('\t'.join([s['id'], s['external_id']]) + '\n')
                utils.gsutil_cp(sample_map_local_fpath, sample_map_bucket_path)
                annotate_job = dataproc.hail_dataproc_job(
                    b,
                    f'batch_seqr_loader/scripts/make_annotated_mt.py '
                    f'--source-path {found_jc_vcf_path} '
                    f'--dest-mt-path {annotated_mt_path} '
                    f'--bucket {anno_tmp_bucket} '
                    '--disable-validation '
                    '--make-checkpoints '
                    f'--remap-tsv {sample_map_bucket_path} '
                    + (f'--vep-block-size {vep_block_size} ' if vep_block_size else ''),
                    max_age='16h',
                    packages=utils.DATAPROC_PACKAGES,
                    num_secondary_workers=utils.NUMBER_OF_DATAPROC_WORKERS,
                    job_name=f'Annotate {project}',
                    vep='GRCh38',
                    depends_on=[jc_job] if jc_job else [],
                )

        if end_with_stage == 'annotate':
            logger.info(
                f'Latest stage is {end_with_stage}, not creating ES index for {project}'
            )
            continue

        timestamp = time.strftime("%Y%m%d-%H%M%S")
        dataproc.hail_dataproc_job(
            b,
            f'batch_seqr_loader/scripts/load_to_es.py '
            f'--mt-path {annotated_mt_path} '
            f'--es-index {project}-{output_version}-{timestamp} '
            f'--es-index-min-num-shards 1 '
            f'--genome-version GRCh38 '
            f'{"--prod" if prod else ""}',
            max_age='16h',
            packages=utils.DATAPROC_PACKAGES,
            num_secondary_workers=10,
            job_name=f'{project}: add to the ES index',
            depends_on=[annotate_job],
            scopes=['cloud-platform'],
        )


if __name__ == '__main__':
    main()  # pylint: disable=E1120
