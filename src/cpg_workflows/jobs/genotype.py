"""
CRAM to GVCF: create Hail Batch jobs to genotype individual sequencing groups.
"""

import logging

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config, image_path, reference_path
from cpg_utils.hail_batch import command, fasta_res_group
from cpg_workflows import utils
from cpg_workflows.filetypes import CramPath, GvcfPath
from cpg_workflows.resources import HIGHMEM, STANDARD

from .picard import get_intervals


def genotype(
    b: hb.Batch,
    sequencing_group_name: str,
    tmp_prefix: Path,
    cram_path: CramPath,
    job_attrs: dict[str, str] | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> list[Job]:
    """
    Takes a CRAM file and runs GATK tools to make a GVCF.
    """
    hc_gvcf_path = tmp_prefix / 'haplotypecaller' / f'{sequencing_group_name}.g.vcf.gz'
    if utils.can_reuse(output_path, overwrite):
        return []

    # collect all the real jobs - [reuse] jobs substituted for None
    jobs = [
        job
        for job in haplotype_caller(
            b=b,
            sequencing_group_name=sequencing_group_name,
            job_attrs=job_attrs,
            output_path=hc_gvcf_path,
            cram_path=cram_path,
            tmp_prefix=tmp_prefix,
            scatter_count=get_config()['workflow'].get('scatter_count_genotype', 50),
            overwrite=overwrite,
            dragen_mode=dragen_mode,
        )
        if job is not None
    ]
    postproc_j = postproc_gvcf(
        b=b,
        gvcf_path=GvcfPath(hc_gvcf_path),
        sequencing_group_name=sequencing_group_name,
        job_attrs=job_attrs,
        output_path=output_path,
        overwrite=overwrite,
    )

    # only keep elements which are not None
    # if both exist, set the dependency
    if postproc_j and jobs:
        postproc_j.depends_on(*jobs)

    if postproc_j:
        jobs.append(postproc_j)

    return jobs


intervals: list[hb.ResourceFile] | None = None


def haplotype_caller(
    b: hb.Batch,
    sequencing_group_name: str,
    cram_path: CramPath,
    tmp_prefix: Path,
    scatter_count: int,
    job_attrs: dict[str, str] | None = None,
    output_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> list[Job | None]:
    """
    Run GATK Haplotype Caller in parallel, split by intervals.
    """
    if utils.can_reuse(output_path, overwrite):
        logging.info(f'Reusing HaplotypeCaller {output_path}, job attrs: {job_attrs}')
        return [None]

    jobs: list[Job | None] = []

    if scatter_count > 1:
        global intervals
        if intervals is None:
            intervals_j, intervals = get_intervals(
                b=b,
                scatter_count=scatter_count,
                source_intervals_path=get_config()['workflow'].get('intervals_path'),
                job_attrs=job_attrs,
                output_prefix=tmp_prefix / f'intervals_{scatter_count}',
            )
            if intervals_j:
                jobs.append(intervals_j)

        hc_fragments = []
        # Splitting variant calling by intervals
        for idx in range(scatter_count):
            assert intervals[idx], intervals
            # give each fragment a tmp location
            fragment = tmp_prefix / 'haplotypecaller' / f'{idx}_of_{scatter_count}_{sequencing_group_name}.g.vcf.gz'
            j, result = _haplotype_caller_one(
                b,
                sequencing_group_name=sequencing_group_name,
                cram_path=cram_path,
                job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
                interval=intervals[idx],
                out_gvcf_path=fragment,
                dragen_mode=dragen_mode,
                overwrite=overwrite,
            )
            hc_fragments.append(result)

            # only consider jobs which weren't scheduled for [reuse]
            if j is not None:
                jobs.append(j)

        merge_j = merge_gvcfs_job(
            b=b,
            sequencing_group_name=sequencing_group_name,
            gvcf_groups=hc_fragments,
            job_attrs=job_attrs,
            out_gvcf_path=output_path,
            overwrite=overwrite,
        )
        if merge_j:
            jobs.append(merge_j)
    else:
        hc_j, _result = _haplotype_caller_one(
            b,
            sequencing_group_name=sequencing_group_name,
            job_attrs=job_attrs,
            cram_path=cram_path,
            out_gvcf_path=output_path,
            overwrite=overwrite,
        )
        if hc_j:
            jobs.append(hc_j)

    return jobs


def _haplotype_caller_one(
    b: hb.Batch,
    sequencing_group_name: str,
    cram_path: CramPath,
    job_attrs: dict | None = None,
    interval: hb.Resource | None = None,
    out_gvcf_path: Path | None = None,
    overwrite: bool = False,
    dragen_mode: bool = True,
) -> tuple[Job, hb.ResourceGroup] | tuple[None, str]:
    """
    Add one GATK HaplotypeCaller job on an interval.
    """
    job_name = 'HaplotypeCaller'

    if utils.can_reuse(out_gvcf_path, overwrite):
        logging.info(f'Reusing HaplotypeCaller {out_gvcf_path}, job attrs: {job_attrs}')
        return None, str(out_gvcf_path)

    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='gatk HaplotypeCaller'))
    j.image(image_path('gatk'))

    # Enough storage to localize CRAMs (can't pass GCS URL to CRAM to gatk directly
    # because we will hit GCP egress bandwidth limit:
    # https://batch.hail.populationgenomics.org.au/batches/7493/jobs/2)
    # CRAMs can be as big as 80G:
    # https://batch.hail.populationgenomics.org.au/batches/74042/jobs/3346
    # plus we need enough space to fit output GVCF (1G) and reference data (5G).
    # HaplotypeCaller is not parallelized, so we request the minimal possible chunk
    # of a Hail worker (2 cores), plus with the `highmem` instance that would
    # give enough memory: 13G. That's not going to give enough disk storage, so we
    # are explicitly requesting more storage.
    #
    # Based on an audit of RD crams on 19/05/23, 99% of crams are <34Gb. Will set the
    # default to 40Gb for genomes then use a run specific confg to run the rare
    # sequencing group that will fail from this limit.
    if (haplo_storage := get_config()['resource_overrides'].get('haplotypecaller_storage')) is not None:
        storage_gb = haplo_storage
    elif get_config()['workflow']['sequencing_type'] == 'genome':
        storage_gb = 40
    else:
        storage_gb = None  # avoid extra disk for exomes

    job_res = HIGHMEM.request_resources(ncpu=2)
    # enough for input CRAM and output GVCF
    job_res.attach_disk_storage_gb = storage_gb
    job_res.set_to_job(j)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sequencing_group_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sequencing_group_name + '.g.vcf.gz.tbi',
        },
    )

    reference = fasta_res_group(b)

    assert isinstance(j.output_gvcf, hb.ResourceGroup)

    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{sequencing_group_name}.cram
    CRAI=$BATCH_TMPDIR/{sequencing_group_name}.cram.crai

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    gatk --java-options \
    "{job_res.java_mem_options()} \
    -XX:GCTimeLimit=50 \
    -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R {reference.base} \\
    -I $CRAM \\
    --read-index $CRAI \\
    {f"-L {interval} " if interval is not None else ""} \\
    --disable-spanning-event-genotyping \\
    {"--dragen-mode " if dragen_mode else ""} \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
    -ERC GVCF \\
    --create-output-variant-index
    """
    j.command(command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j, j.output_gvcf


def merge_gvcfs_job(
    b: hb.Batch,
    sequencing_group_name: str,
    gvcf_groups: list[str | hb.ResourceGroup],
    job_attrs: dict | None = None,
    out_gvcf_path: Path | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Combine by-interval GVCFs into a single sequencing group wide GVCF file.
    """
    job_name = f'Merge {len(gvcf_groups)} GVCFs'
    if utils.can_reuse(out_gvcf_path, overwrite):
        logging.info(f'Reusing {job_name} {out_gvcf_path}, job attrs: {job_attrs}')
        return None

    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='picard MergeVcfs'))

    j.image(image_path('picard'))
    j.cpu(2)
    java_mem = 11
    j.memory('highmem')  # ~ 6G/core ~ 12G
    j.storage(f'{len(gvcf_groups) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sequencing_group_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sequencing_group_name + '.g.vcf.gz.tbi',
        },
    )

    input_cmd = ''
    for gvcf_group in gvcf_groups:
        # if the output was recoverable, read into the batch
        if isinstance(gvcf_group, str):
            gvcf_group = b.read_input_group(**{'g.vcf.gz': gvcf_group, 'g.vcf.gz.tbi': f'{gvcf_group}.tbi'})
        input_cmd += f'INPUT={gvcf_group["g.vcf.gz"]} '

    assert isinstance(j.output_gvcf, hb.ResourceGroup)
    cmd = f"""
    picard -Xms{java_mem}g \
    MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}
    """
    j.command(command(cmd))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


def postproc_gvcf(
    b: hb.Batch,
    gvcf_path: GvcfPath,
    sequencing_group_name: str,
    job_attrs: dict | None = None,
    overwrite: bool = False,
    output_path: Path | None = None,
    depends_on: list[Job] | None = None,
) -> Job | None:
    """
    1. Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration, and reduce the number of GVCF blocking bins to 2.
    2. Subsets GVCF to main, not-alt chromosomes to avoid downstream errors.
    3. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    4. Renames the GVCF sequencing group name to use CPG ID.
    """
    logging.info(f'Adding GVCF postproc job for sequencing group {sequencing_group_name}, gvcf {gvcf_path}')
    job_name = 'Postproc GVCF'
    if utils.can_reuse(output_path, overwrite):
        logging.info(f'Reusing {output_path} output for {job_name}. {job_attrs}')
        return None

    j = b.new_job(job_name, (job_attrs or {}) | dict(tool='gatk ReblockGVCF'))
    j.image(image_path('gatk'))

    # We need at least 2 CPU, so on 16-core instance it would be 8 jobs,
    # meaning we have more than enough disk (265/8=33.125G).
    # Enough to fit a pre-reblocked GVCF, which can be as big as 10G,
    # the reblocked result (1G), and ref data (5G).
    storage_gb = get_config()['resource_overrides'].get('postproc_gvcf_storage', 20)
    job_res = STANDARD.set_resources(j, ncpu=2, storage_gb=storage_gb)

    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        },
    )

    reference = fasta_res_group(b)
    noalt_regions = b.read_input(reference_path('broad/noalt_bed'))
    gvcf = b.read_input(str(gvcf_path.path))
    gq_bands = get_config()['workflow']['reblock_gq_bands']

    assert isinstance(j.output_gvcf, hb.ResourceGroup)

    cmd = f"""\
    GVCF={gvcf}
    GVCF_NODP=$BATCH_TMPDIR/{sequencing_group_name}-nodp.g.vcf.gz
    REBLOCKED=$BATCH_TMPDIR/{sequencing_group_name}-reblocked.g.vcf.gz

    # Reindexing just to make sure the index is not corrupted
    bcftools index --tbi $GVCF

    # Remove INFO/DP field, which contradicts the FORMAT/DP, in the way that
    # it has _all_ reads, not just variant-calling-usable reads. If we keep INFO/DP,
    # ReblockGVCF would prioritize it over FORMAT/DP and change FORMAT/DP to INFO/DP
    # in the resulting merged blocks. It would pick the highest INFO/DP when merging
    # multiple blocks, so a variant in a small homopolymer region (surrounded by
    # long DP=0 areas), that attracted piles of low-MQ reads with INFO/DP=1000
    # will translate into a long GQ<20 block with the same FORMAT/DP=1000,
    # which is wrong, because most of this block has no reads.
    bcftools view $GVCF \\
    | bcftools annotate -x INFO/DP \\
    | bcftools view -Oz -o $GVCF_NODP
    tabix -p vcf $GVCF_NODP

    gatk --java-options "{job_res.java_mem_options()}" \\
    ReblockGVCF \\
    --reference {reference.base} \\
    -V $GVCF_NODP \\
    -do-qual-approx \\
    --floor-blocks {' '.join(f'-GQB {b}' for b in gq_bands)} \\
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools query -l $GVCF)

    bcftools view $REBLOCKED -T {noalt_regions} \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "$EXISTING_SN {sequencing_group_name}") \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    tabix -p vcf {j.output_gvcf['g.vcf.gz']}
    """
    j.command(command(cmd, setup_gcp=True, monitor_space=True, define_retry_function=True))
    if output_path:
        b.write_output(j.output_gvcf, str(output_path).replace('.g.vcf.gz', ''))
    if depends_on:
        j.depends_on(*depends_on)
    return j
