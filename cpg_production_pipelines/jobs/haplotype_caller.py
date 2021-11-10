import logging
from os.path import join
from typing import Optional, List

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_production_pipelines import resources, utils
from cpg_production_pipelines.jobs import wrap_command, new_job
from cpg_production_pipelines.pipeline import Batch
from cpg_production_pipelines.smdb import SMDB
from cpg_production_pipelines.jobs import split_intervals

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def produce_gvcf(
    b: hb.Batch,
    output_path: str,
    sample_name: str,
    project_name: str,
    tmp_bucket: str,
    cram_path: str,
    crai_path: Optional[str] = None,
    number_of_intervals: Optional[int] = 1,
    intervals: Optional[hb.ResourceGroup] = None,
    overwrite: bool = True,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    external_id: Optional[str] = None,
    dragen_mode: bool = False,
) -> Job:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    j = new_job(b, 'make GVCF', sample_name, project_name)
    if utils.can_reuse(output_path, overwrite):
        j.name += ' [reuse]'
        return j

    logger.info(
        f'Submitting the variant calling jobs to write {output_path} for {sample_name}'
    )

    reference = b.read_input_group(
        base=resources.REF_FASTA,
        fai=resources.REF_FASTA + '.fai',
        dict=resources.REF_FASTA.replace('.fasta', '')
        .replace('.fna', '')
        .replace('.fa', '')
        + '.dict',
    )

    hc_gvcf_path = join(tmp_bucket, 'haplotypecaller', f'{sample_name}.g.vcf.gz')
    hc_jobs = []
    if number_of_intervals is not None and number_of_intervals > 1:
        if intervals is None:
            intervals = split_intervals.get_intervals(
                b=b,
                scatter_count=number_of_intervals,
                out_bucket=join(tmp_bucket, 'intervals'),
            )

        # Splitting variant calling by intervals
        for idx in range(number_of_intervals):
            hc_jobs.append(
                hc_job(
                    b,
                    sample_name=sample_name,
                    project_name=project_name,
                    reference=reference,
                    cram_fpath=cram_path,
                    crai_fpath=crai_path,
                    interval=intervals[f'interval_{idx}'],
                    interval_idx=idx,
                    number_of_intervals=number_of_intervals,
                    depends_on=depends_on,
                    dragen_mode=dragen_mode,
                )
            )
        hc_j = merge_gvcfs_job(
            b=b,
            sample_name=sample_name,
            project_name=project_name,
            gvcfs=[j.output_gvcf for j in hc_jobs],
            out_gvcf_path=hc_gvcf_path,
            overwrite=overwrite,
        )
    else:
        hc_j = hc_job(
            b,
            sample_name=sample_name,
            project_name=project_name,
            reference=reference,
            cram_fpath=cram_path,
            crai_fpath=crai_path,
            depends_on=depends_on,
            out_gvcf_path=hc_gvcf_path,
            overwrite=overwrite,
        )
        hc_jobs.append(hc_j)

    if depends_on:
        hc_jobs[0].depends_on(*depends_on)

    postproc_j = postproc_gvcf(
        b=b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf_path=hc_gvcf_path,
        out_gvcf_path=output_path,
        reference=reference,
        overwrite=overwrite,
        depends_on=[hc_j],
        external_id=external_id,
    )
    last_j = postproc_j

    if smdb:
        last_j = smdb.add_running_and_completed_update_jobs(
            b=b,
            analysis_type='gvcf',
            output_path=output_path,
            sample_names=[sample_name],
            project_name=project_name,
            first_j=hc_jobs[0],
            last_j=last_j,
            depends_on=depends_on,
        )
    return last_j


def hc_job(
    b: Batch,
    sample_name: str,
    project_name: str,
    reference: hb.ResourceGroup,
    cram_fpath: str,
    crai_fpath: Optional[str] = None,
    interval: Optional[hb.ResourceFile] = None,
    interval_idx: Optional[int] = None,
    number_of_intervals: int = 1,
    depends_on: Optional[List[Job]] = None,
    out_gvcf_path: Optional[str] = None,
    overwrite: bool = True,
    dragen_mode: bool = False,
) -> Job:
    """
    Run HaplotypeCaller on an input BAM or CRAM, and output GVCF
    """
    job_name = 'HaplotypeCaller'
    if interval_idx is not None:
        job_name += f', {interval_idx}/{number_of_intervals}'

    j = new_job(b, job_name, sample_name, project_name)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.GATK_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage('60G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    if depends_on:
        j.depends_on(*depends_on)

    cmd = f"""\
    gatk --java-options "-Xms{java_mem}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R {reference.base} \\
    -I {cram_fpath} \\
    --read-index {crai_fpath or (cram_fpath + '.crai')} \\
    {f"-L {interval} " if interval is not None else ""} \\
    --disable-spanning-event-genotyping \\
    {"--dragen-mode " if dragen_mode else ""} \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 20 \\
    -ERC GVCF
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def merge_gvcfs_job(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    gvcfs: List[hb.ResourceGroup],
    out_gvcf_path: Optional[str],
    overwrite: bool = True,
) -> Job:
    """
    Combine by-interval GVCFs into a single sample GVCF file
    """
    job_name = f'Merge {len(gvcfs)} GVCFs'
    j = new_job(b, job_name, sample_name, project_name)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j
    
    j.image(resources.SAMTOOLS_PICARD_IMAGE)
    j.cpu(2)
    java_mem = 7
    j.memory('standard')  # ~ 4G/core ~ 7.5G
    j.storage(f'{len(gvcfs) * 1.5 + 2}G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    input_cmd = ' '.join(f'INPUT={g["g.vcf.gz"]}' for g in gvcfs)
    cmd = f"""
    picard -Xms{java_mem}g \
    MergeVcfs {input_cmd} OUTPUT={j.output_gvcf['g.vcf.gz']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def postproc_gvcf(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf_path: str,
    reference: hb.ResourceGroup,
    out_gvcf_path: str,
    overwrite: bool,
    depends_on: Optional[List[Job]] = None,
    external_id: Optional[str] = None,
) -> Job:
    if utils.can_reuse(out_gvcf_path, overwrite):
        return new_job(b, 'Postproc GVCF [reuse]', sample_name, project_name)

    logger.info(
        f'Adding reblock and subset jobs for sample {sample_name}, gvcf {out_gvcf_path}'
    )
    reblock_j = reblock_gvcf(
        b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf=b.read_input_group(
            **{'g.vcf.gz': input_gvcf_path, 'g.vcf.gz.tbi': input_gvcf_path + '.tbi'}
        ),
        reference=reference,
        overwrite=overwrite,
    )
    if depends_on:
        reblock_j.depends_on(*depends_on)
    subset_to_noalt_j = subset_noalt(
        b,
        sample_name=sample_name,
        project_name=project_name,
        input_gvcf=reblock_j.output_gvcf,
        noalt_regions=b.read_input(resources.NOALT_REGIONS),
        overwrite=overwrite,
        out_gvcf_path=out_gvcf_path,
        external_sample_id=external_id,
        internal_sample_id=sample_name,
    )
    return subset_to_noalt_j


def reblock_gvcf(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    overwrite: bool = True,
    out_gvcf_path: Optional[str] = None,
) -> Job:
    """
    Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration
    """
    j = new_job(b, 'ReblockGVCF', sample_name, project_name)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    cmd = f"""\
    gatk --java-options "-Xms{mem_gb - 1}g" ReblockGVCF \\
    --reference {reference.base} \\
    -V {input_gvcf['g.vcf.gz']} \\
    -do-qual-approx \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    --create-output-variant-index true
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j


def subset_noalt(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    input_gvcf: hb.ResourceGroup,
    noalt_regions: str,
    overwrite: bool,
    out_gvcf_path: Optional[str] = None,
    external_sample_id: Optional[str] = None,
    internal_sample_id: Optional[str] = None,
) -> Job:
    """
    1. Subset GVCF to main chromosomes to avoid downstream errors
    2. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    3. Renames sample name from external_sample_id to internal_sample_id
    """
    j = new_job(b, 'SubsetToNoalt', sample_name, project_name)
    if utils.can_reuse(out_gvcf_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(resources.BCFTOOLS_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'30G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )
    if external_sample_id and internal_sample_id:
        reheader_cmd = f"""
        | bcftools reheader -s <(echo "{external_sample_id} {internal_sample_id}")
        """
    else:
        reheader_cmd = ''

    cmd = f"""
    bcftools view {input_gvcf['g.vcf.gz']} -T {noalt_regions} \\
    | bcftools annotate -x INFO/DS \\
    {reheader_cmd} \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    bcftools index --tbi {j.output_gvcf['g.vcf.gz']}
    """
    j.command(wrap_command(cmd, monitor_space=True))
    if out_gvcf_path:
        b.write_output(j.output_gvcf, out_gvcf_path.replace('.g.vcf.gz', ''))
    return j
