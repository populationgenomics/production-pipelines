"""
Create Hail Batch jobs for seqr_loader SNPs and Indels cohort annotation.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import command, get_batch, query_command
from cpg_workflows.query_modules import seqr_loader
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse

from .picard import get_intervals


def annotate_cohort_jobs_snps_indels(
    vcf_path: Path,
    out_mt_path: Path,
    vep_ht_path: Path,
    checkpoint_prefix: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Annotate cohort VCF for seqr loader, SNPs and Indels.
    """

    j = get_batch().new_job('Annotate cohort', job_attrs)
    j.image(config_retrieve['workflow', 'driver_image'])
    j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(vcf_path),
            str(out_mt_path),
            str(vep_ht_path),
            str(checkpoint_prefix),
        ),
    )
    return j


def split_vcf_for_vep(
    b: hb.Batch,
    merged_vcf_path: Path,
    tmp_bucket: Path,
    out_siteonly_vcf_path: Path,
    out_siteonly_vcf_part_paths: list[Path] | None = None,
    intervals_path: Path | None = None,
    exclude_intervals_path: Path | None = None,
    job_attrs: dict | None = None,
    scatter_count: int | None = None,
) -> list[Job]:
    """
    Takes the merged VCF from the seqr_loader_long_read and prepares it
    for annotation with VEP by splitting it into parts and creating site-only VCFs.
    """
    jobs: list[Job] = []
    intervals: list[str] | list[hb.ResourceFile] = []
    if intervals_path and intervals_path.suffix == '.bed':
        # If intervals_path is a bed file, read the intervals directly
        intervals_j = None
        intervals = get_intervals_from_bed(intervals_path)
        assert scatter_count == len(intervals)
    else:
        # If intervals_path is not specified, use the get_intervals picard job
        intervals_j, intervals = get_intervals(
            b=b,
            source_intervals_path=intervals_path,
            exclude_intervals_path=exclude_intervals_path,
            scatter_count=scatter_count,
            job_attrs=job_attrs,
            output_prefix=tmp_bucket / f'intervals_{scatter_count}',
        )
    if intervals_j:
        jobs.append(intervals_j)

    all_output_paths = [out_siteonly_vcf_path]
    if out_siteonly_vcf_part_paths:
        assert len(out_siteonly_vcf_part_paths) == scatter_count
        all_output_paths.extend(out_siteonly_vcf_part_paths)
    if can_reuse(all_output_paths + [to_path(f'{p}.tbi') for p in all_output_paths]):
        return []

    vcfs: list[hb.ResourceGroup] = []
    siteonly_vcfs: list[hb.ResourceGroup] = []
    for idx, interval in enumerate(intervals):
        out_vcf_path = tmp_bucket / 'split-merged-vcf' / 'parts' / f'part{idx + 1}.vcf.gz'
        if out_siteonly_vcf_part_paths:
            siteonly_jc_vcf_path = out_siteonly_vcf_part_paths[idx]
        else:
            siteonly_jc_vcf_path = (
                (tmp_bucket / 'split-merged-vcf-siteonly' / 'parts' / f'part{idx + 1}.vcf.gz')
                if scatter_count > 1
                else out_siteonly_vcf_path
            )
        vcf_j, vcf = _add_split_vcf_job(
            b=b,
            input_vcf=merged_vcf_path,
            output_vcf_path=out_vcf_path,
            interval=interval,
            job_attrs=job_attrs,
        )
        if vcf_j:
            jobs.append(vcf_j)
            vcfs.append(vcf)

        siteonly_j, siteonly_j_vcf = add_make_sitesonly_job(
            b=b,
            input_vcf=vcfs[idx],
            output_vcf_path=siteonly_jc_vcf_path,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        # make add_make_sitesonly_job return the vcf file as apointer
        siteonly_vcfs.append(siteonly_j_vcf)
        if siteonly_j:
            jobs.append(siteonly_j)

    jobs = [j for j in jobs if j is not None]
    for j in jobs:
        j.name = f'Long Read SNPs Indels Annotation: {j.name}'
    return jobs


def _add_split_vcf_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    interval: str | hb.ResourceFile,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
) -> hb.ResourceGroup:
    """
    Split VCF by interval and write to part VCF.
    Use GATK SelectVariants to split the VCF by interval.
    """
    if can_reuse(output_vcf_path):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            },
        )
    job_name = 'SplitVcf'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk SelectVariants'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2)
    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    assert isinstance(j.output_vcf, hb.ResourceGroup)

    cmd = f"""\
    gatk --java-options "{res.java_mem_options()}" \\
    SelectVariants \\
    -V {input_vcf['vcf.gz']} \\
    -O {j.output_vcf['vcf.gz']} \\
    -L {interval}
    """
    j.command(command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True))
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def add_make_sitesonly_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
    storage_gb: int | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Create sites-only VCF with only site-level annotations.

    Returns: a Job object and a single output ResourceGroup j.sites_only_vcf
    """
    if output_vcf_path and can_reuse(output_vcf_path):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            },
        )

    job_name = 'MakeSitesOnlyVcf'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk MakeSitesOnlyVcf'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2)
    if storage_gb:
        j.storage(f'{storage_gb}G')
    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    assert isinstance(j.output_vcf, hb.ResourceGroup)
    j.command(
        command(
            f"""
    gatk --java-options "{res.java_mem_options()}" \\
    MakeSitesOnlyVcf \\
    -I {input_vcf['vcf.gz']} \\
    -O {j.output_vcf['vcf.gz']}

    if [[ ! -e {j.output_vcf['vcf.gz.tbi']} ]]; then
        tabix -p vcf {j.output_vcf['vcf.gz']}
    fi
    """,
        ),
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def get_intervals_from_bed(intervals_path: Path) -> list[str]:
    """
    Read intervals from a bed file.
    Increment the start position of each interval by 1 to match the 1-based
    coordinate system used by GATK.
    """
    with intervals_path.open('r') as f:
        intervals = []
        for line in f:
            chrom, start, end = line.strip().split('\t')
            intervals.append(f'{chrom}:{int(start)+1}-{end}')
    return intervals
