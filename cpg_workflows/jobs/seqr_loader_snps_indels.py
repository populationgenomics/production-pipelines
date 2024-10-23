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
from cpg_workflows.utils import can_reuse, get_intervals_from_bed

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
    j.image(config_retrieve(['workflow', 'driver_image']))
    j.command(
        query_command(
            seqr_loader,
            seqr_loader.annotate_cohort.__name__,
            str(vcf_path),
            str(out_mt_path),
            str(vep_ht_path),
            None,  # site_only_vqsr_vcf_path
            str(checkpoint_prefix),
            True,  # long_read
        ),
    )
    return j


def split_merged_vcf_and_get_sitesonly_vcfs_for_vep(
    b: hb.Batch,
    scatter_count: int,
    merged_vcf_path: Path,
    tmp_bucket: Path,
    out_siteonly_vcf_part_paths: list[Path] | None = None,
    intervals_path: Path | None = None,
    exclude_intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Takes the merged VCF from the seqr_loader_long_read pipeline and prepares it
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

    # all_output_paths = [out_siteonly_vcf_path]
    all_output_paths = []
    if out_siteonly_vcf_part_paths:
        assert len(out_siteonly_vcf_part_paths) == scatter_count
        all_output_paths.extend(out_siteonly_vcf_part_paths)
    if can_reuse(all_output_paths + [to_path(f'{p}.tbi') for p in all_output_paths]):
        return []

    merged_vcf = b.read_input_group(
        **{
            'vcf.gz': str(merged_vcf_path),
            'vcf.gz.tbi': str(merged_vcf_path) + '.tbi',
        },
    )

    split_vcfs_paths = [
        tmp_bucket / 'split-merged-vcf' / 'parts' / f'part{idx + 1}.vcf.gz' for idx in range(scatter_count)
    ]
    siteonly_vcfs: list[hb.ResourceGroup] = []

    split_vcf_j = add_split_vcf_job(
        b=b,
        input_vcf=merged_vcf,
        intervals=intervals,
        output_vcf_paths=split_vcfs_paths,
        job_attrs=(job_attrs or {}) | dict(part='all'),
    )

    for idx in range(scatter_count):
        if out_siteonly_vcf_part_paths:
            siteonly_vcf_path = out_siteonly_vcf_part_paths[idx]
        else:
            siteonly_vcf_path = tmp_bucket / 'split-merged-vcf-siteonly' / 'parts' / f'part{idx + 1}.vcf.gz'
        vcf_part = b.read_input_group(
            **{
                'vcf.gz': str(split_vcfs_paths[idx]),
                'vcf.gz.tbi': str(split_vcfs_paths[idx]) + '.tbi',
            },
        )
        siteonly_j, siteonly_j_vcf = add_make_sitesonly_job(
            b=b,
            input_vcf=vcf_part,
            output_vcf_path=siteonly_vcf_path,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )

        siteonly_vcfs.append(siteonly_j_vcf)
        if siteonly_j:
            jobs.append(siteonly_j)
            if split_vcf_j:
                siteonly_j.depends_on(split_vcf_j)

    jobs = [j for j in jobs if j is not None]
    for j in jobs:
        j.name = f'Long Read SNPs Indels Annotation: {j.name}'
    return jobs


def add_split_vcf_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    intervals: list[str] | list[hb.ResourceFile],
    output_vcf_paths: list[Path],
    job_attrs: dict | None = None,
) -> Job:
    """
    Split VCF by interval and write to part VCF.
    Use GATK SelectVariants to split the VCF by interval.
    """
    j = b.new_job('SplitVcf', (job_attrs or {}) | {'tool': 'gatk SelectVariants'})
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2)

    for idx, interval in enumerate(intervals):
        output_vcf_path = output_vcf_paths[idx]
        j.declare_resource_group(
            **{
                str(idx): {
                    'vcf.gz': '{root}.vcf.gz',
                    'vcf.gz.tbi': '{root}.vcf.gz.tbi',
                },
            },
        )
        cmd = f"""\
            gatk --java-options "{res.java_mem_options()}" \\
            SelectVariants \\
            -V {input_vcf['vcf.gz']} \\
            -O {j[str(idx)]['vcf.gz']} \\
            -L {interval}
        """
        j.command(command(cmd, monitor_space=False, setup_gcp=True, define_retry_function=True))
        b.write_output(j[str(idx)], str(output_vcf_path).replace('.vcf.gz', ''))

    # Wait for all parts to be written before returning
    j.command('wait && echo "All parts written"')

    return j


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
