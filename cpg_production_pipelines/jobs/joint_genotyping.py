import json
import logging
from enum import Enum
from os.path import join, basename
from typing import Optional, List, Collection, Dict, Tuple, Set

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines import resources, utils
from cpg_production_pipelines.jobs import wrap_command
from cpg_production_pipelines.jobs import split_intervals
from cpg_production_pipelines.pipeline import Sample
from cpg_production_pipelines.smdb import SMDB

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class JointGenotyperTool(Enum):
    GenotypeGVCFs = 1
    GnarlyGenotyper = 2


def make_joint_genotyping_jobs(
    b: hb.Batch,
    out_jc_vcf_path: str,
    samples: Collection[Sample],
    genomicsdb_bucket: str,
    tmp_bucket: str,
    gvcf_by_sid: Dict[str, str],
    local_tmp_dir: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    tool: JointGenotyperTool = JointGenotyperTool.GnarlyGenotyper,
    do_filter_excesshet: bool = True
) -> Job:
    """
    Assumes all samples have a 'file' of 'type'='gvcf' in `samples_df`.
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    job_name = 'Joint-calling+VQSR'
    if utils.can_reuse([out_jc_vcf_path], overwrite):
        return b.new_job(f'{job_name} [reuse]')

    logger.info(f'Submitting the joint-calling and VQSR jobs.')

    # For small callsets, we don't apply the ExcessHet filtering anyway
    do_filter_excesshet = len(samples) >= 1000 and do_filter_excesshet

    genomicsdb_path_per_interval = dict()
    for idx in range(resources.NUMBER_OF_GENOMICS_DB_INTERVALS):
        genomicsdb_path_per_interval[idx] = join(
            genomicsdb_bucket,
            f'interval_{idx}_outof_{resources.NUMBER_OF_GENOMICS_DB_INTERVALS}',
        )
    # Determining which samples to add. Using the first interval, so the assumption
    # is that all DBs have the same set of samples.
    (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        sample_names_to_remove,
        updating_existing_db,
        sample_map_bucket_path,
    ) = _samples_to_add_to_db(
        genomicsdb_gcs_path=genomicsdb_path_per_interval[0],
        interval_idx=0,
        local_tmp_dir=local_tmp_dir,
        samples=samples,
        tmp_bucket=tmp_bucket,
        gvcf_by_sid=gvcf_by_sid,
    )
    sample_ids = set(s.id for s in samples)
    assert sample_names_will_be_in_db == sample_ids
    samples_hash = utils.hash_sample_ids(sample_ids)

    intervals = split_intervals.get_intervals(
        b=b,
        scatter_count=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
    )

    import_gvcfs_job_per_interval = dict()
    first_j = None
    if sample_names_to_add:
        logger.info(f'Queueing genomics-db-import jobs')
        for idx in range(resources.NUMBER_OF_GENOMICS_DB_INTERVALS):
            import_gvcfs_job, _ = _add_import_gvcfs_job(
                b=b,
                genomicsdb_gcs_path=genomicsdb_path_per_interval[idx],
                sample_names_to_add=sample_names_to_add,
                sample_names_to_skip=sample_names_already_added,
                sample_names_to_remove=sample_names_to_remove,
                sample_names_will_be_in_db=sample_names_will_be_in_db,
                updating_existing_db=updating_existing_db,
                sample_map_bucket_path=sample_map_bucket_path,
                interval=intervals[f'interval_{idx}'],
                interval_idx=idx,
                number_of_intervals=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
            )
            import_gvcfs_job_per_interval[idx] = import_gvcfs_job
            if first_j is None:
                first_j = import_gvcfs_job

    scattered_vcf_by_interval: Dict[int, hb.ResourceGroup] = dict()

    joint_calling_tmp_bucket = f'{tmp_bucket}/joint_calling/{samples_hash}'
    for idx in range(resources.NUMBER_OF_GENOMICS_DB_INTERVALS):
        joint_called_vcf_path = (
            f'{joint_calling_tmp_bucket}/by_interval/interval_{idx}.vcf.gz'
        )
        if utils.can_reuse(joint_called_vcf_path, overwrite):
            j = b.new_job('Joint genotyping [reuse]')
            scattered_vcf_by_interval[idx] = b.read_input_group(
                **{
                    'vcf.gz': joint_called_vcf_path,
                    'vcf.gz.tbi': joint_called_vcf_path + '.tbi',
                }
            )
            if first_j is None:
                first_j = j
        else:
            logger.info(f'Queueing {tool.name} job')
            genotype_vcf_job = _add_joint_genotyper_job(
                b,
                genomicsdb_path=genomicsdb_path_per_interval[idx],
                overwrite=overwrite,
                number_of_samples=len(sample_names_will_be_in_db),
                interval_idx=idx,
                number_of_intervals=resources.NUMBER_OF_GENOMICS_DB_INTERVALS,
                interval=intervals[f'interval_{idx}'],
                tool=tool,
            )
            if first_j is None:
                first_j = genotype_vcf_job
            if import_gvcfs_job_per_interval.get(idx):
                genotype_vcf_job.depends_on(import_gvcfs_job_per_interval.get(idx))

            if do_filter_excesshet:
                logger.info(f'Queueing exccess het filter job')
                exccess_filter_job = _add_exccess_het_filter(
                    b,
                    input_vcf=genotype_vcf_job.output_vcf,
                    interval=intervals[f'interval_{idx}'],
                    overwrite=overwrite,
                )
                last_job = exccess_filter_job
            else:
                last_job = genotype_vcf_job
            
            last_job = _add_make_sites_only_job(
                b=b,
                input_vcf=last_job.output_vcf,
                overwrite=overwrite,
            )
            scattered_vcf_by_interval[idx] = last_job.output_vcf

    scattered_vcfs = list(scattered_vcf_by_interval.values())
    logger.info(f'Queueing gather-VCF job')
    final_gathered_vcf_job = _add_final_gather_vcf_job(
        b,
        input_vcfs=scattered_vcfs,
        overwrite=overwrite,
        output_vcf_path=out_jc_vcf_path,
        is_allele_specific=True,
    )
    if first_j is None:
        first_j = final_gathered_vcf_job
    last_j = final_gathered_vcf_job

    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type='joint-calling',
            output_path=out_jc_vcf_path,
            sample_names=sample_ids,
            first_j=first_j,
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def _samples_to_add_to_db(
    genomicsdb_gcs_path,
    interval_idx,
    local_tmp_dir,
    samples,
    tmp_bucket: str,
    gvcf_by_sid: Dict[str, str],
) -> Tuple[Set[str], Set[str], Set[str], Set[str], bool, str]:
    if utils.file_exists(join(genomicsdb_gcs_path, 'callset.json')):
        # Checking if samples exists in the DB already
        genomicsdb_metadata = join(local_tmp_dir, f'callset-{interval_idx}.json')
        utils.gsutil_cp(
            src_path=join(genomicsdb_gcs_path, 'callset.json'),
            dst_path=genomicsdb_metadata,
            disable_check_hashes=True,
        )
        with open(genomicsdb_metadata) as f:
            db_metadata = json.load(f)
        sample_names_in_db = set(s['sample_name'] for s in db_metadata['callsets'])
        sample_names_requested = set([s.id for s in samples])
        sample_names_to_add = sample_names_requested - sample_names_in_db
        sample_names_to_remove = sample_names_in_db - sample_names_requested
        if sample_names_to_remove:
            # GenomicsDB doesn't support removing, so creating a new DB
            updating_existing_db = False
            sample_names_already_added = set()
            sample_names_to_add = {s.id for s in samples}
            sample_names_will_be_in_db = sample_names_to_add
            logger.info(
                f'GenomicDB {genomicsdb_gcs_path} exists, but '
                f'{len(sample_names_to_remove)} samples need '
                f'to be removed: {", ".join(sample_names_to_remove)}, so creating a new '
                f'DB with {len(sample_names_will_be_in_db)} samples: '
                f'{", ".join(sample_names_will_be_in_db)}'
            )
        else:
            updating_existing_db = True
            sample_names_will_be_in_db = sample_names_in_db | sample_names_to_add
            sample_names_already_added = sample_names_requested & sample_names_in_db
            if sample_names_already_added:
                logger.info(
                    f'{len(sample_names_already_added)} samples '
                    f'{", ".join(sample_names_already_added)} already exist in the DB '
                    f'{genomicsdb_gcs_path}, skipping adding them.'
                )
            if sample_names_to_remove:
                logger.info(
                    f'There are {len(sample_names_to_remove)} samples that need to be '
                    f'removed from the DB {genomicsdb_gcs_path}: '
                    f'{", ".join(sample_names_to_remove)}. Re-creating the DB '
                    f'using the updated set of samples'
                )
            elif sample_names_to_add:
                logger.info(
                    f'Will add {len(sample_names_to_add)} samples '
                    f'{", ".join(sample_names_to_add)} into the DB {genomicsdb_gcs_path}'
                )
            else:
                logger.warning(
                    f'Nothing will be added into the DB {genomicsdb_gcs_path}'
                )
    else:
        # Initiate new DB
        sample_names_already_added = set()
        sample_names_to_add = {s.id for s in samples}
        sample_names_will_be_in_db = sample_names_to_add
        sample_names_to_remove = set()
        updating_existing_db = False
        logger.info(
            f'GenomicDB {genomicsdb_gcs_path} doesn\'t exist, so creating a new one '
            f'with {len(sample_names_to_add)} samples: {", ".join(sample_names_to_add)}'
        )

    sample_map_bucket_path = join(tmp_bucket, 'genomicsdb', 'sample_map.csv')
    sample_map_local_fpath = join(local_tmp_dir, basename(sample_map_bucket_path))
    with open(sample_map_local_fpath, 'w') as f:
        for sid in sample_names_to_add:
            f.write('\t'.join([sid, gvcf_by_sid[sid]]) + '\n')
    utils.gsutil_cp(sample_map_local_fpath, sample_map_bucket_path)

    return (
        sample_names_to_add,
        sample_names_will_be_in_db,
        sample_names_already_added,
        sample_names_to_remove,
        updating_existing_db,
        sample_map_bucket_path,
    )


def _add_import_gvcfs_job(
    b: hb.Batch,
    genomicsdb_gcs_path: str,
    sample_names_to_add: Set[str],
    sample_names_to_skip: Set[str],
    sample_names_to_remove: Set[str],
    sample_names_will_be_in_db: Set[str],
    updating_existing_db: bool,
    sample_map_bucket_path: str,
    interval: hb.ResourceFile,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
) -> Tuple[Optional[Job], Set[str]]:
    """
    Add GVCFs to a genomics database (or create a new instance if it doesn't exist)
    Returns a Job, or None if no new samples to add
    """
    rm_cmd = ''

    if updating_existing_db:
        # Update existing DB
        genomicsdb_param = f'--genomicsdb-update-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Adding to GenomicsDB'
        msg = f'Adding {len(sample_names_to_add)} samples: {", ".join(sample_names_to_add)}'
        if sample_names_to_skip:
            msg += (
                f'Skipping adding {len(sample_names_to_skip)} samples that are already '
                f'in the DB: {", ".join(sample_names_to_skip)}"'
            )
    else:
        # Initiate new DB
        genomicsdb_param = f'--genomicsdb-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Creating GenomicsDB'
        if sample_names_to_remove:
            # Need to remove the existing database
            rm_cmd = f'gsutil -q rm -rf {genomicsdb_gcs_path}'
        msg = (
            f'Creating a new DB with {len(sample_names_to_add)} samples: '
            f'{", ".join(sample_names_to_add)}'
        )

    sample_map = b.read_input(sample_map_bucket_path)

    if interval_idx is not None:
        job_name += f' {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(resources.GATK_IMAGE)
    ncpu = 16
    j.cpu(ncpu)
    java_mem = 16
    j.memory('lowmem')  # ~ 1G/core ~ 14.4G
    if depends_on:
        j.depends_on(*depends_on)

    j.command(wrap_command(f"""\
    # We've seen some GenomicsDB performance regressions related to intervals, 
    # so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg. There's no data in between since 
    # we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    
    # The batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!

    echo "{msg}"

    {rm_cmd}

    gatk --java-options -Xms{java_mem}g \\
      GenomicsDBImport \\
      {genomicsdb_param} \\
      --batch-size 50 \\
      -L {interval} \\
      --sample-name-map {sample_map} \\
      --reader-threads {ncpu} \\
      --merge-input-intervals \\
      --consolidate

    """, monitor_space=True, setup_gcp=True))
    return j, sample_names_will_be_in_db


def _add_joint_genotyper_job(
    b: hb.Batch,
    genomicsdb_path: str,
    overwrite: bool,
    number_of_samples: int,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
    tool: JointGenotyperTool = JointGenotyperTool.GnarlyGenotyper,
) -> Job:
    """
    Runs GATK GnarlyGenotyper or GenotypeGVCFs on a combined_gvcf VCF bgzipped file,
    pre-called with HaplotypeCaller.

    GenotypeGVCFs is a standard GATK joint-genotyping tool.

    GnarlyGenotyper is an experimental GATK joint-genotyping tool that performs 
    "quick and dirty" joint genotyping on large cohorts, of GVCFs post-processed 
    with ReblockGVCF.

    HaplotypeCaller must be used with `-ERC GVCF` or `-ERC BP_RESOLUTION` to add
    genotype likelihoods.

    ReblockGVCF must be run to add all the annotations necessary for VQSR:
    QUALapprox, VarDP, RAW_MQandDP.
    """
    job_name = f'Joint genotyping: {tool.name}'
    if interval_idx is not None:
        job_name += f' {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name)
    if utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    reference = b.read_input_group(
        base=resources.REF_FASTA,
        fai=resources.REF_FASTA + '.fai',
        dict=resources.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
    )

    j.image(resources.GATK_IMAGE)
    j.cpu(2)
    j.memory('standard')  # ~ 4G/core ~ 8G
    # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals
    j.storage(f'{4 + number_of_samples * 4 // number_of_intervals}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    cmd = f"""\
    gatk --java-options -Xms8g \\
    {tool.name} \\
    -R {reference.base} \\
    -O {j.output_vcf['vcf.gz']} \\
    -D {resources.DBSNP_VCF} \\
    -V gendb.{genomicsdb_path} \\
    {f'-L {interval} ' if interval else ''} \\
    --only-output-calls-starting-in-intervals \\
    """
    if tool == JointGenotyperTool.GnarlyGenotyper:
        cmd += f"""\
    --keep-all-sites \\
    --create-output-variant-index
    """
    else:
        cmd += f"""\
    --merge-input-intervals \\
    -G AS_StandardAnnotation
    """
    j.command(wrap_command(cmd, monitor_space=True, setup_gcp=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j


def _add_exccess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int = 32,
    excess_het_threshold: float = 54.69,
    interval: Optional[hb.ResourceGroup] = None,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Filter a large cohort callset on Excess Heterozygosity.

    The filter applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinuity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: ExcessHet filter'
    j = b.new_job(job_name)
    if utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(wrap_command(f"""
    # Captring stderr to avoid Batch pod from crashing with OOM from millions of
    # warning messages from VariantFiltration, e.g.:
    # > JexlEngine - ![0,9]: 'ExcessHet > 54.69;' undefined variable ExcessHet
    gatk --java-options -Xms3g \\
    VariantFiltration \\
    --filter-expression 'ExcessHet > {excess_het_threshold}' \\
    --filter-name ExcessHet \\
    {f'-L {interval} ' if interval else ''} \\
    -O {j.output_vcf['vcf.gz']} \\
    -V {input_vcf['vcf.gz']} \\
    2> {j.stderr}
    """
    ))
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_make_sites_only_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int = 32,
    output_vcf_path: Optional[str] = None,
) -> Job:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    j = b.new_job(job_name)
    if utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.GATK_IMAGE)
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    j.command(wrap_command(f"""
    gatk --java-options -Xms6g \\
    MakeSitesOnlyVcf \\
    -I {input_vcf['vcf.gz']} \\
    -O {j.output_vcf['vcf.gz']}
    """))
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))

    return j


def _add_final_gather_vcf_job(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    overwrite: bool,
    output_vcf_path: str = None,
    is_allele_specific: bool = True,
) -> Job:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    job_name = f'Gather {len(input_vcfs)} VCFs'
    j = b.new_job(job_name)
    if utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{1 + len(input_vcfs) * (0.1 if is_allele_specific else 2)}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(wrap_command(f"""
    # --ignore-safety-checks makes a big performance difference so we include it in 
    # our invocation. This argument disables expensive checks that the file headers 
    # contain the same set of genotyped samples and that files are in order 
    # by position of first record.
    gatk --java-options -Xms{java_mem}g \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}

    tabix {j.output_vcf['vcf.gz']}
    """, monitor_space=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, output_vcf_path.replace('.vcf.gz', ''))
    return j
