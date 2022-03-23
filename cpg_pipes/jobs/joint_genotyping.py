"""
Create Hail Batch jobs for joint genotyping.
"""

import json
import logging
from enum import Enum

import pandas as pd
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import Path
from cpg_pipes import images, utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs import split_intervals
from cpg_pipes.jobs.vcf import gather_vcfs
from cpg_pipes.pipeline.targets import Sample
from cpg_pipes.types import GvcfPath
from cpg_pipes.refdata import RefData

logger = logging.getLogger(__file__)


class JointGenotyperTool(Enum):
    """
    Tool used for joint genotyping. GenotypeGVCFs is more stable, 
    GnarlyGenotyper is fater but more experimental.
    """
    GenotypeGVCFs = 1
    GnarlyGenotyper = 2
    

def make_joint_genotyping_jobs(
    b: hb.Batch,
    out_vcf_path: Path,
    out_siteonly_vcf_path: Path,
    samples: list[Sample],
    genomicsdb_bucket: Path,
    tmp_bucket: Path,
    refs: RefData,
    gvcf_by_sid: dict[str, GvcfPath],
    overwrite: bool = True,
    # Default to GenotypeGVCFs because Gnarly is a bit weird, e.g. it adds <NON_REF>
    # variants with AC_adj annotations (other variants have AC):
    # bcftools view gs://cpg-fewgenomes-test/unittest/inputs/chr20/gnarly/joint-called-siteonly.vcf.gz | zgrep 7105364
    tool: JointGenotyperTool = JointGenotyperTool.GenotypeGVCFs,
    do_filter_excesshet: bool = True,
    scatter_count: int | None = None,
) -> list[Job]:
    """
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    job_name = 'Joint genotyping'
    if utils.can_reuse([out_vcf_path, out_siteonly_vcf_path], overwrite):
        return [b.new_job(f'{job_name} [reuse]')]
   
    if len(samples) == 0:
        raise ValueError(
            'Provided samples collection for joint calling should contain '
            'at least one active sample'
        )

    sequencing_types = set(s.sequencing_type for s in samples)
    if len(sequencing_types) > 1:
        raise ValueError(
            f'Samples of more than one sequencing type used for joint calling: '
            f'{sequencing_types}. Joint calling can be only performed on samples '
            f'of the same type.'
        )
    sequencing_type = list(sequencing_types)[0]

    logger.info(f'Submitting joint-calling jobs.')

    scatter_count = scatter_count or RefData.number_of_joint_calling_intervals

    jobs: list[Job] = []
    
    intervals_j = split_intervals.get_intervals(
        b=b,
        refs=refs,
        sequencing_type=sequencing_type,
        scatter_count=scatter_count,
    )
    intervals = [intervals_j[f'intervals{i}.list'] for i in range(scatter_count)]
    jobs.append(intervals_j)
    
    # There are some problems with using GenomcsDB in cloud (`genomicsdb_cloud()`):
    # Using cloud + --consolidate on larger datasets causes a TileDB error:
    # https://github.com/broadinstitute/gatk/issues/7653
    # Also, disabling --consolidate decreases performance of reading, causing
    # GenotypeGVCFs to run for too long. So for now it's safer to disable cloud,
    # and use the version of function that passes a tarball around (`genomicsdb()`):
    import_gvcfs_job_per_interval, genomicsdb_path_per_interval = genomicsdb(
        b=b,
        samples=samples,
        genomicsdb_bucket=genomicsdb_bucket,
        tmp_bucket=tmp_bucket,
        gvcf_by_sid=gvcf_by_sid,
        intervals=intervals,
        scatter_count=scatter_count,
        overwrite=overwrite,
    )
    jobs.extend([v for k, v in import_gvcfs_job_per_interval.items()])

    vcf_by_interval: dict[int, hb.ResourceGroup] = dict()
    siteonly_vcf_by_interval: dict[int, hb.ResourceGroup] = dict()

    sample_ids = set(s.id for s in samples)
    samples_hash = utils.hash_sample_ids(list(sample_ids))
    jc_tmp_bucket = tmp_bucket / 'joint_calling' / samples_hash
    for idx in range(scatter_count):
        jc_vcf_path = jc_tmp_bucket / 'by_interval' / f'interval_{idx}.vcf.gz'
        filt_jc_vcf_path = (
            jc_tmp_bucket / 'by_interval_excess_het_filter' / f'interval_{idx}.vcf.gz'
        )
        siteonly_jc_vcf_path = (
           jc_tmp_bucket / 'by_interval_site_only' / f'interval_{idx}.vcf.gz'
        )

        jc_vcf_j, jc_vcf = _add_joint_genotyper_job(
            b,
            genomicsdb_path=genomicsdb_path_per_interval[idx],
            overwrite=overwrite,
            number_of_samples=len(samples),
            refs=refs,
            interval_idx=idx,
            number_of_intervals=scatter_count,
            interval=intervals[idx],
            tool=tool,
            output_vcf_path=jc_vcf_path,
        )
        vcf_by_interval[idx] = jc_vcf
        jobs.append(jc_vcf_j)
        if import_gvcfs_job_per_interval.get(idx):
            jc_vcf_j.depends_on(import_gvcfs_job_per_interval.get(idx))

        # For small callsets, we don't apply the ExcessHet filtering anyway
        if len(samples) >= 1000 and do_filter_excesshet:
            logger.info(f'Queueing exccess het filter job')
            exccess_filter_j, exccess_filter_jc_vcf = _add_exccess_het_filter(
                b,
                input_vcf=jc_vcf,
                interval=intervals[idx],
                overwrite=overwrite,
                output_vcf_path=filt_jc_vcf_path,
            )
            vcf_by_interval[idx] = exccess_filter_jc_vcf
            jobs.append(exccess_filter_j)

        siteonly_j, siteonly_vcf = _add_make_sitesonly_job(
            b=b,
            input_vcf=vcf_by_interval[idx],
            overwrite=overwrite,
            output_vcf_path=siteonly_jc_vcf_path,
        )
        siteonly_vcf_by_interval[idx] = siteonly_vcf
        jobs.append(siteonly_j)

    logger.info(f'Queueing gather VCFs job')
    gather_j, _ = gather_vcfs(
        b,
        input_vcfs=[v for k, v in vcf_by_interval.items()],
        overwrite=overwrite,
        output_vcf_path=out_vcf_path,
        site_only=False,
    )
    gather_j.name = 'Joint genotyping: ' + gather_j.name
    jobs.append(gather_j)

    logger.info(f'Queueing gather site-only VCFs job')
    gather_siteonly_j, _ = gather_vcfs(
        b,
        input_vcfs=[v for k, v in siteonly_vcf_by_interval.items()],
        overwrite=overwrite,
        output_vcf_path=out_siteonly_vcf_path,
        site_only=True,
    )
    gather_siteonly_j.name = 'Joint genotyping: ' + gather_siteonly_j.name
    jobs.append(gather_siteonly_j)
    return jobs


def genomicsdb(
    b: hb.Batch,
    samples: list[Sample],
    genomicsdb_bucket: Path,
    tmp_bucket: Path,
    gvcf_by_sid: dict[str, GvcfPath],
    intervals: list[hb.Resource],
    scatter_count: int = RefData.number_of_joint_calling_intervals,    
    overwrite: bool = False,
) -> tuple[dict[int, Job], dict[int, Path]]:
    """
    Create GenomicDBs for each interval, given new samples.
    """
    sample_map_bucket_path = tmp_bucket / 'genomicsdb' / 'sample_map.csv'
    df = pd.DataFrame([{
        'id': s.id, 
        'path': str(gvcf_by_sid[s.id].path)
    } for s in samples])
    df.to_csv(str(sample_map_bucket_path), index=False, header=False, sep='\t')

    path_per_interval = {
        idx: genomicsdb_bucket / f'interval_{idx}_outof_{scatter_count}.tar' 
        for idx in range(scatter_count)
    }

    job_per_interval = dict()
    for idx in range(scatter_count):
        out_path = path_per_interval[idx]

        job_name = 'Joint genotyping: creating GenomicsDB'
        if idx is not None:
            job_name += f' {idx + 1}/{scatter_count}'
        j = b.new_job(job_name)
        job_per_interval[idx] = j

        if utils.can_reuse(out_path, overwrite):
            j.name += ' [reuse]'
            continue

        j.image(images.GATK_IMAGE)
    
        sample_map = b.read_input(str(sample_map_bucket_path))
        
        # The Broad: testing has shown that the multithreaded reader initialization
        # does not scale well beyond 5 threads, so don't increase beyond that.
        nthreads = 5
    
        # The Broad: The memory setting here is very important and must be several
        # GiB lower than the total memory allocated to the VM because this tool uses
        # a significant amount of non-heap memory for native libraries.
        xms_gb = 8
        xmx_gb = 25
    
        STANDARD.set_resources(
            j,
            nthreads=nthreads, 
            mem_gb=xmx_gb + 1,
            storage_gb=20,
        )

        params = [
            # The Broad: We've seen some GenomicsDB performance regressions related 
            # to intervals, so we're going to pretend we only have a single interval
            # using the --merge-input-intervals arg. There's no data in between since we
            # didn't run HaplotypeCaller over those loci, so we're not wasting any compute
            '--merge-input-intervals',
            '--consolidate',
            # The batch_size value was carefully chosen here as it is the optimal value for 
            # the amount of memory allocated within the task; please do not change it 
            # without consulting the Hellbender (GATK engine) team!
            '--batch-size', '50',
        ]

        cmd = f"""\
        WORKSPACE={out_path.with_suffix('').name}
        
        gatk --java-options "-Xms{xms_gb}g -Xmx{xmx_gb}g" \\
        GenomicsDBImport \\
        --genomicsdb-workspace-path $WORKSPACE \\
        -L {intervals[idx]} \\
        --sample-name-map {sample_map} \\
        --reader-threads {nthreads} \\
        {" ".join(params)}
    
        tar -cf {j.db_tar} $WORKSPACE
        """
        j.command(wrap_command(cmd, monitor_space=True, setup_gcp=True))
        b.write_output(j.db_tar, str(out_path))
    return job_per_interval, path_per_interval


def genomicsdb_cloud(
    b: hb.Batch,
    samples: list[Sample],
    genomicsdb_bucket: Path,
    tmp_bucket: Path,
    gvcf_by_sid: dict[str, GvcfPath],
    intervals: hb.ResourceGroup,
    scatter_count: int = RefData.number_of_joint_calling_intervals,    
    depends_on: list[Job] | None = None,
) -> tuple[dict[int, Job], dict[int, Path]]:
    """
    Update or create GenomicDBs for each interval, given new samples.
    """
    genomicsdb_path_per_interval = dict()
    for idx in range(scatter_count):
        genomicsdb_path_per_interval[idx] = (
            genomicsdb_bucket / f'interval_{idx}_outof_{scatter_count}'
        )
        
    # Determining which samples to add. Using the first interval, so the assumption
    # is that all DBs have the same set of samples.
    (
        sample_names_to_add,
        sample_names_already_added,
        sample_names_to_remove,
        updating_existing_db,
        sample_map_bucket_path,
    ) = _samples_to_add_to_db(
        genomicsdb_gcs_path=genomicsdb_path_per_interval[0],
        samples=samples,
        tmp_bucket=tmp_bucket,
        gvcf_by_sid=gvcf_by_sid,
    )

    import_gvcfs_job_per_interval = dict()
    if sample_names_to_add:
        logger.info(f'Queueing genomics-db-import jobs')
        for idx in range(scatter_count):
            import_gvcfs_job, _ = _genomicsdb_import_cloud(
                b=b,
                genomicsdb_gcs_path=genomicsdb_path_per_interval[idx],
                sample_names_to_add=sample_names_to_add,
                sample_names_to_skip=sample_names_already_added,
                sample_names_will_be_in_db=set(s.id for s in samples),
                updating_existing_db=updating_existing_db,
                sample_map_bucket_path=sample_map_bucket_path,
                interval=intervals[f'interval_{idx}'],
                interval_idx=idx,
                number_of_intervals=scatter_count,
            )
            if depends_on:
                import_gvcfs_job.depends_on(*depends_on)
            import_gvcfs_job_per_interval[idx] = import_gvcfs_job
    return import_gvcfs_job_per_interval, genomicsdb_path_per_interval


def _samples_to_add_to_db(
    genomicsdb_gcs_path: Path,
    samples: list[Sample],
    tmp_bucket: Path,
    gvcf_by_sid: dict[str, GvcfPath],
) -> tuple[set[str], set[str], set[str], bool, Path]:
    if utils.exists(genomicsdb_gcs_path / 'callset.json'):
        # Checking if samples exists in the DB already
        with (genomicsdb_gcs_path / 'callset.json').open() as f:
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
            assert sample_names_will_be_in_db == sample_names_requested
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

    sample_map_bucket_path = tmp_bucket / 'genomicsdb' / 'sample_map.csv'
    with sample_map_bucket_path.open('w') as f:
        for sid in sample_names_to_add:
            f.write('\t'.join([sid, str(gvcf_by_sid[sid].path)]) + '\n')

    assert sample_names_will_be_in_db == {s.id for s in samples}
    return (
        sample_names_to_add,
        sample_names_already_added,
        sample_names_to_remove,
        updating_existing_db,
        sample_map_bucket_path,
    )


def _genomicsdb_import_cloud(
    b: hb.Batch,
    genomicsdb_gcs_path: Path,
    sample_names_to_add: set[str],
    sample_names_to_skip: set[str],
    sample_names_will_be_in_db: set[str],
    updating_existing_db: bool,
    sample_map_bucket_path: Path,
    interval: hb.ResourceFile,
    interval_idx: int | None = None,
    number_of_intervals: int = 1,
    depends_on: list[Job] | None = None,
    overwrite: bool = False,
) -> tuple[Job, set[str]]:
    """
    Add GVCFs to a genomics database (or create a new instance if it doesn't exist)
    Returns a Job, or None if no new samples to add
    """
    rm_cmd = ''
    
    if updating_existing_db and not overwrite:
        # Update existing DB
        genomicsdb_param = f'--genomicsdb-update-workspace-path {genomicsdb_gcs_path}'
        job_name = 'Joint genotyping: adding to GenomicsDB'
        msg = f'Adding {len(sample_names_to_add)} samples: {", ".join(sample_names_to_add)}'
        if sample_names_to_skip:
            msg += (
                f'Skipping adding {len(sample_names_to_skip)} samples that are already '
                f'in the DB: {", ".join(sample_names_to_skip)}"'
            )
    else:
        # Initiate new DB
        job_name = 'Joint genotyping: creating GenomicsDB'

        genomicsdb_param = f'--genomicsdb-workspace-path {genomicsdb_gcs_path}'
        # Need to remove the existing database if exists
        rm_cmd = f'gsutil ls {genomicsdb_gcs_path} && gsutil -q rm -rf {genomicsdb_gcs_path}'
        msg = (
            f'Creating a new DB with {len(sample_names_to_add)} samples: '
            f'{", ".join(sample_names_to_add)}'
        )

    sample_map = b.read_input(str(sample_map_bucket_path))

    if interval_idx is not None:
        job_name += f' {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name)
    j.image(images.GATK_IMAGE)

    if depends_on:
        j.depends_on(*depends_on)
    
    # The Broad: testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    nthreads = 5

    # The Broad: The memory setting here is very important and must be several
    # GiB lower than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    xms_gb = 8
    xmx_gb = 25

    STANDARD.set_resources(j, nthreads=nthreads, mem_gb=xmx_gb + 1)
    
    params = [
        # The Broad: We've seen some GenomicsDB performance regressions related 
        # to intervals, so we're going to pretend we only have a single interval
        # using the --merge-input-intervals arg. There's no data in between since we
        # didn't run HaplotypeCaller over those loci, so we're not wasting any compute
        '--merge-input-intervals',
        '--consolidate',
        # The batch_size value was carefully chosen here as it is the optimal value for 
        # the amount of memory allocated within the task; please do not change it 
        # without consulting the Hellbender (GATK engine) team!
        '--batch-size', '50',
        # Using cloud + --consolidate causes a TileDB error:
        # https://github.com/broadinstitute/gatk/issues/7653
        # However disabling --consolidate decreases performance of reading, cause
        # GenotypeGVCFs to run forever. So instead, using non-cloud GenomicsDB.
    ]

    j.command(wrap_command(f"""\

    echo "{msg}"

    {rm_cmd}

    gatk --java-options "-Xms{xms_gb}g -Xmx{xmx_gb}g" \
    GenomicsDBImport \\
    {genomicsdb_param} \\
    -L {interval} \\
    --sample-name-map {sample_map} \\
    --reader-threads {nthreads} \\
    {" ".join(params)}
    """, monitor_space=True, setup_gcp=True))
    return j, sample_names_will_be_in_db


def _add_joint_genotyper_job(
    b: hb.Batch,
    genomicsdb_path: Path,
    overwrite: bool,
    number_of_samples: int,
    refs: RefData,
    interval_idx: int | None = None,
    number_of_intervals: int = 1,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    tool: JointGenotyperTool = JointGenotyperTool.GnarlyGenotyper,
) -> tuple[Job, hb.ResourceGroup]:
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
        return j, b.read_input_group(**{
            'vcf.gz': str(output_vcf_path),
            'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
        })

    j.image(images.GATK_IMAGE)

    xms_gb = 8
    xmx_gb = 25

    STANDARD.set_resources(
        j, 
        mem_gb=xmx_gb + 1,
        # 4G (fasta+fai+dict) + 4G per sample divided by the number of intervals:
        storage_gb=4 + number_of_samples * 4 // number_of_intervals
    )

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    reference = refs.fasta_res_group(b)
    
    if str(genomicsdb_path).endswith('.tar'):
        # can't use directly from cloud, need to copy and uncompress:
        input_cmd = f"""\
        retry_gs_cp {genomicsdb_path} {genomicsdb_path.name}
        tar -xf {genomicsdb_path.name}
        WORKSPACE=gendb://{genomicsdb_path.with_suffix('').name}
        """
    else:
        input_cmd = f"""\
        WORKSPACE=gendb.{genomicsdb_path}
        """

    cmd = f"""\
    {input_cmd}
    
    gatk --java-options "-Xms{xms_gb}g -Xmx{xmx_gb}g" \\
    {tool.name} \\
    -R {reference.base} \\
    -O {j.output_vcf['vcf.gz']} \\
    -D {refs.dbsnp_vcf} \\
    -V $WORKSPACE \\
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
    j.command(wrap_command(
        cmd, monitor_space=True, setup_gcp=True, 
        define_retry_function=True
    ))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def _add_exccess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int = 32,
    excess_het_threshold: float = 54.69,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
) -> tuple[Job, hb.ResourceGroup]:
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
        return j, b.read_input_group(**{
            'vcf.gz': str(output_vcf_path),
            'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
        })

    j.image(images.GATK_IMAGE)
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
    """))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def _add_make_sitesonly_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int = 32,
    output_vcf_path: Path | None = None,
) -> tuple[Job, hb.ResourceGroup]:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    j = b.new_job(job_name)
    if output_vcf_path and utils.can_reuse(output_vcf_path, overwrite):
        j.name += ' [reuse]'
        return j, b.read_input_group(**{
            'vcf.gz': str(output_vcf_path),
            'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
        })

    j.image(images.GATK_IMAGE)
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
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf
