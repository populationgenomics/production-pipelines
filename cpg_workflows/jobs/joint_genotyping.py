"""
Create Hail Batch jobs for joint genotyping.
"""

import logging
import math
from enum import Enum

import pandas as pd
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path, command
from cpg_utils import Path
import hailtop.batch as hb
from hailtop.batch import Resource
from hailtop.batch.job import Job

from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import GvcfPath

from .vcf import gather_vcfs
from .picard import get_intervals


class JointGenotyperTool(Enum):
    """
    Tool used for joint genotyping. GenotypeGVCFs is more stable,
    GnarlyGenotyper is faster but more experimental.
    """

    GenotypeGVCFs = 1
    GnarlyGenotyper = 2


def make_joint_genotyping_jobs(
    b: hb.Batch,
    gvcf_by_sid: dict[str, GvcfPath],
    out_vcf_path: Path,
    out_siteonly_vcf_path: Path,
    tmp_bucket: Path,
    overwrite: bool = False,
    # Default to GenotypeGVCFs because Gnarly is a bit weird, e.g. it adds <NON_REF>
    # variants with AC_adj annotations (other variants have AC):
    # bcftools view gs://cpg-fewgenomes-test/unittest/inputs/chr20/gnarly/joint-called-siteonly.vcf.gz | zgrep 7105364
    tool: JointGenotyperTool = JointGenotyperTool.GenotypeGVCFs,
    do_filter_excesshet: bool = True,
    intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Adds samples to the GenomicsDB and runs joint genotyping on them.
    Outputs a multi-sample VCF under `output_vcf_path`.
    """
    if can_reuse([out_vcf_path, out_siteonly_vcf_path], overwrite):
        return []
    if len(gvcf_by_sid) == 0:
        raise ValueError(
            'Provided samples collection for joint calling should contain '
            'at least one active sample'
        )
    logging.info(f'Submitting joint-calling jobs')

    scatter_count = 50
    if len(gvcf_by_sid) > 300:
        scatter_count = 100
    if len(gvcf_by_sid) > 1000:
        scatter_count = 200

    jobs: list[Job] = []
    intervals_j, intervals = get_intervals(
        b=b,
        source_intervals_path=intervals_path,
        scatter_count=scatter_count,
        job_attrs=job_attrs,
        output_prefix=tmp_bucket / f'intervals_{scatter_count}',
    )
    if intervals_j:
        jobs.append(intervals_j)

    vcfs: list[hb.ResourceFile] = []
    siteonly_vcfs: list[hb.ResourceFile] = []

    # Preparing inputs for GenomicsDB
    genomicsdb_bucket = tmp_bucket / 'genomicsdbs'
    sample_map_bucket_path = genomicsdb_bucket / 'sample_map.csv'
    df = pd.DataFrame(
        [{'id': sid, 'path': str(path)} for sid, path in gvcf_by_sid.items()]
    )
    if not get_config()['workflow'].get('dry_run', False):
        with sample_map_bucket_path.open('w') as fp:
            df.to_csv(fp, index=False, header=False, sep='\t')

    do_filter_excesshet = len(gvcf_by_sid) >= 1000 and do_filter_excesshet

    for idx, interval in enumerate(intervals):
        genomicsdb_path = (
            genomicsdb_bucket / f'interval_{idx + 1}_outof_{scatter_count}.tar'
        )
        import_gvcfs_j = genomicsdb(
            b=b,
            sample_map_bucket_path=sample_map_bucket_path,
            interval=interval,
            output_path=genomicsdb_path,
            overwrite=overwrite,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )

        jc_vcf_path = (
            (tmp_bucket / 'joint-genotyper' / 'parts' / f'part{idx + 1}.vcf.gz')
            if scatter_count > 1
            else out_vcf_path
        )
        filt_jc_vcf_path = (
            (tmp_bucket / 'excess-filter' / 'parts' / f'part{idx + 1}.vcf.gz')
            if scatter_count > 1 and do_filter_excesshet
            else out_vcf_path
        )
        siteonly_jc_vcf_path = (
            (tmp_bucket / 'siteonly' / 'parts' / f'part{idx + 1}.vcf.gz')
            if scatter_count > 1
            else out_siteonly_vcf_path
        )

        jc_vcf_j, jc_vcf = _add_joint_genotyper_job(
            b,
            genomicsdb_path=genomicsdb_path,
            overwrite=overwrite,
            number_of_samples=len(gvcf_by_sid),
            interval=interval,
            tool=tool,
            output_vcf_path=jc_vcf_path,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        if jc_vcf_j:
            jobs.append(jc_vcf_j)
            if import_gvcfs_j:
                jc_vcf_j.depends_on(import_gvcfs_j)

        # For small callsets, we don't apply the ExcessHet filtering anyway
        if len(gvcf_by_sid) >= 1000 and do_filter_excesshet:
            logging.info(f'Queueing excess het filter job')
            excess_filter_j, excess_filter_jc_vcf = _add_excess_het_filter(
                b,
                input_vcf=jc_vcf,
                interval=interval,
                overwrite=overwrite,
                output_vcf_path=filt_jc_vcf_path,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            )
            if excess_filter_j:
                jobs.append(excess_filter_j)
                if jc_vcf_j:
                    excess_filter_j.depends_on(jc_vcf_j)
            vcfs.append(excess_filter_jc_vcf['vcf.gz'])
        else:
            vcfs.append(jc_vcf['vcf.gz'])

        siteonly_j, siteonly_vcf = _add_make_sitesonly_job(
            b=b,
            input_vcf=vcfs[idx],
            overwrite=overwrite,
            output_vcf_path=siteonly_jc_vcf_path,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        siteonly_vcfs.append(siteonly_vcf['vcf.gz'])
        if siteonly_j:
            jobs.append(siteonly_j)

    if scatter_count > 1:
        logging.info(f'Queueing gather VCFs job')
        gather_j, _ = gather_vcfs(
            b,
            input_vcfs=vcfs,
            overwrite=overwrite,
            out_vcf_path=out_vcf_path,
            site_only=False,
            gvcf_count=len(gvcf_by_sid),
            job_attrs=job_attrs,
        )
        if gather_j:
            gather_j.name = f'Joint genotyping: {gather_j.name}'
            jobs.append(gather_j)

        logging.info(f'Queueing gather site-only VCFs job')
        gather_siteonly_j, _ = gather_vcfs(
            b,
            input_vcfs=siteonly_vcfs,
            overwrite=overwrite,
            out_vcf_path=out_siteonly_vcf_path,
            site_only=True,
            job_attrs=job_attrs,
        )
        if gather_siteonly_j:
            gather_siteonly_j.name = f'Joint genotyping: {gather_siteonly_j.name}'
            jobs.append(gather_siteonly_j)

    return jobs


def genomicsdb(
    b: hb.Batch,
    sample_map_bucket_path: Path,
    output_path: Path,
    interval: Resource | None = None,
    overwrite: bool = False,
    job_attrs: dict | None = None,
) -> Job | None:
    """
    Create GenomicDBs for an interval.

    Note that GenomicsDBImport and GenotypeGVCFs can interact directly with a DB
    sitting on a bucket, without having to make a copy. However, there some problems
    with that setup: writing/reading DB from a bucket together with the `--consolidate`
    option on larger datasets would cause this TileDB error:
    https://github.com/broadinstitute/gatk/issues/7653
    And disabling `--consolidate` would decrease the performance of reading, causing
    GenotypeGVCFs to run for too long, making it much more expensive. So for now, it's
    safer to copy the DB between instances in a tarball. In the future, we might adopt
    Hail Query VCF combiner for all this.
    """
    if can_reuse(output_path, overwrite):
        return None

    job_name = 'Joint genotyping: creating GenomicsDB'
    j = b.new_job(
        job_name,
        (job_attrs or {})
        | dict(
            tool='gatk GenomicsDBImport',
        ),
    )
    j.image(image_path('gatk'))

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
        # The Broad:
        # > We've seen some GenomicsDB performance regressions related
        #   to intervals, so we're going to pretend we only have a single interval
        #   using the --merge-input-intervals arg. There's no data in between since we
        #   didn't run HaplotypeCaller over those loci, so we're not wasting any
        #   compute.
        '--merge-input-intervals',
        '--consolidate',
        # The Broad:
        # > The batch_size value was carefully chosen here as it is the optimal value
        #   for the amount of memory allocated within the task; please do not change
        #   it without consulting the Hellbender (GATK engine) team!
        '--batch-size',
        '50',
    ]

    cmd = f"""\
    WORKSPACE={output_path.with_suffix('').name}

    # Multiple GenomicsDBImport read same GVCFs in parallel, which could lead to 
    # some of the jobs failing reading a GVCF, e.g.:
    # https://batch.hail.populationgenomics.org.au/batches/74388/jobs/45 
    # So wrapping the command in a "retry" call, which would attempt it multiple 
    # times after a timeout using.
    # --overwrite-existing-genomicsdb-workspace is to make sure the $WORKSPACE 
    # directory from a previous attempt is not in the way of a new attempt.
    function run {{
    gatk --java-options "-Xms{xms_gb}g -Xmx{xmx_gb}g" \\
    GenomicsDBImport \\
    --genomicsdb-workspace-path $WORKSPACE \\
    {f'-L {interval}' if interval is not None else ''} \\
    --sample-name-map {sample_map} \\
    --reader-threads {nthreads} \\
    --overwrite-existing-genomicsdb-workspace \\
    {" ".join(params)} && \\
    tar -cf {j.db_tar} $WORKSPACE
    }}
    retry run
    """
    j.command(
        command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
    )
    b.write_output(j.db_tar, str(output_path))
    return j


def joint_calling_storage_gb(n_samples: int, sequencing_type: str) -> int:
    """
    Calculate storage for a joint calling job, to fit a joint-called VCF
    and a genomics db. The required storage grows with a number of samples,
    but sublinearly.
    >>> joint_calling_storage_gb(1, 'genome')
    26
    >>> joint_calling_storage_gb(50, 'genome')
    26
    >>> joint_calling_storage_gb(100, 'genome')
    86
    >>> joint_calling_storage_gb(200, 'genome')
    226
    >>> joint_calling_storage_gb(10000, 'genome')
    486
    """
    disk_gb = 26  # Minimal disk (default for a 4-cpu standard machine)
    if n_samples >= 50:
        # for every 10x samples (starting from 50), add {multiplier}G
        multiplier = {
            'genome': 200,
            'exome': 20,
        }[sequencing_type]
        disk_gb += int(multiplier * math.log(n_samples / 50, 10))
    return min(1000, disk_gb)  # requesting 1T max


def _add_joint_genotyper_job(
    b: hb.Batch,
    genomicsdb_path: Path,
    overwrite: bool,
    number_of_samples: int,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    tool: JointGenotyperTool = JointGenotyperTool.GnarlyGenotyper,
    job_attrs: dict | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
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
    `QUALapprox`, `VarDP`, `RAW_MQandDP`.
    """
    if can_reuse(output_vcf_path, overwrite):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            }
        )
    job_name = f'Joint genotyping: {tool.name}'
    job_attrs = (job_attrs or {}) | {'tool': f'gatk {tool.name}'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.request_resources(ncpu=4)
    res.attach_disk_storage_gb = joint_calling_storage_gb(
        number_of_samples,
        get_config()['workflow']['sequencing_type'],
    )
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    reference = fasta_res_group(b)

    if str(genomicsdb_path).endswith('.tar'):
        # can't use directly from cloud, need to copy and uncompress:
        input_cmd = f"""\
        retry_gs_cp {genomicsdb_path} $BATCH_TMPDIR/{genomicsdb_path.name}
        tar -xf $BATCH_TMPDIR/{genomicsdb_path.name} -C $BATCH_TMPDIR/
        WORKSPACE=gendb://$BATCH_TMPDIR/{genomicsdb_path.with_suffix('').name}
        """
    else:
        input_cmd = f"""\
        WORKSPACE=gendb.{genomicsdb_path}
        """

    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""\
    {input_cmd}
    
    gatk --java-options "-Xmx{res.get_java_mem_mb()}m" \\
    {tool.name} \\
    -R {reference.base} \\
    -O {j.output_vcf['vcf.gz']} \\
    -D {reference_path('broad/dbsnp_vcf')} \\
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
    j.command(
        command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def _add_excess_het_filter(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    overwrite: bool,
    disk_size: int = 32,
    excess_het_threshold: float = 54.69,
    interval: hb.Resource | None = None,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Filter a large cohort callset on Excess Heterozygosity.

    The filter applies only to large callsets (`not is_small_callset`)

    Requires all samples to be unrelated.

    ExcessHet estimates the probability of the called samples exhibiting excess
    heterozygosity with respect to the null hypothesis that the samples are unrelated.
    The higher the score, the higher the chance that the variant is a technical artifact
    or that there is consanguinity among the samples. In contrast to Inbreeding
    Coefficient, there is no minimal number of samples for this annotation.

    Returns: a Job object with a single output j.output_vcf of type ResourceGroup
    """
    if can_reuse(output_vcf_path, overwrite):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            }
        )

    job_name = 'Joint genotyping: ExcessHet filter'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    assert isinstance(j.output_vcf, hb.ResourceGroup)
    j.command(
        command(
            f"""
    # Capturing stderr to avoid Batch pod from crashing with OOM from millions of
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
        )
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def _add_make_sitesonly_job(
    b: hb.Batch,
    input_vcf: hb.ResourceFile,
    overwrite: bool,
    disk_size: int = 32,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    if output_vcf_path and can_reuse(output_vcf_path, overwrite):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            }
        )

    job_name = 'Joint genotyping: MakeSitesOnlyVcf'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk MakeSitesOnlyVcf'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    j.memory('8G')
    j.storage(f'{disk_size}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    assert isinstance(j.output_vcf, hb.ResourceGroup)
    j.command(
        command(
            f"""
    gatk --java-options -Xms6g \\
    MakeSitesOnlyVcf \\
    -I {input_vcf} \\
    -O {j.output_vcf['vcf.gz']}
    """
        )
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf
