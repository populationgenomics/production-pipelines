"""
Create Hail Batch jobs for variant calling in individual samples.
"""

import logging
from os.path import join, basename
from typing import Optional, List

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import resources, utils
from cpg_pipes.smdb import SMDB
from cpg_pipes.jobs import split_intervals
from cpg_pipes.hailbatch import wrap_command

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def produce_gvcf(
    b: hb.Batch,
    sample_name: str,
    project_name: str,
    tmp_bucket: str,
    cram_path: str,
    crai_path: Optional[str] = None,
    output_path: Optional[str] = None,
    number_of_intervals: Optional[int] = 1,
    intervals: Optional[hb.ResourceGroup] = None,
    overwrite: bool = True,
    check_existence: bool = True,
    depends_on: Optional[List[Job]] = None,
    smdb: Optional[SMDB] = None,
    dragen_mode: bool = False,
) -> Job:
    """
    Takes all samples with a 'file' of 'type'='bam' in `samples_df`,
    and runs HaplotypeCaller on them, and sets a new 'file' of 'type'='gvcf'

    HaplotypeCaller is run in an interval-based sharded way, with per-interval
    HaplotypeCaller jobs defined in a nested loop.
    """
    if output_path and check_existence and utils.can_reuse(output_path, overwrite):
        return b.new_job('Make GVCF [reuse]', dict(sample=sample_name, project=project_name))

    logger.info(
        f'Submitting the variant calling jobs' + 
        (f' to write {output_path}' if output_path else '') + 
        ' for {sample_name}'
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
                    overwrite=overwrite,
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
        gvcf_path=hc_gvcf_path,
        output_path=output_path,
        overwrite=overwrite,
        depends_on=[hc_j],
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
    b: hb.Batch,
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
        job_name += f', {interval_idx + 1}/{number_of_intervals}'

    j = b.new_job(job_name, dict(sample=sample_name, project=project_name))
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
    j = b.new_job(job_name, dict(sample=sample_name, project=project_name))
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
    gvcf_path: str,
    sample_name: str,
    project_name: str,
    overwrite: bool,
    output_path: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
) -> Job:
    """
    1. Runs ReblockGVCF to annotate with allele-specific VCF INFO fields
    required for recalibration, and reduce the number of GVCF blocking bins to 2.
    2. Subsets GVCF to main, not-alt chromosomes to avoid downstream errors.
    3. Removes the DS INFO field that is added to some HGDP GVCFs to avoid errors
       from Hail about mismatched INFO annotations
    4. Renames the GVCF sample name to use CPG ID.
    """
    if output_path and utils.can_reuse(output_path, overwrite):
        return b.new_job('Postproc GVCF [reuse]', dict(sample=sample_name, project=project_name))

    logger.info(f'Adding GVCF postproc job for sample {sample_name}, gvcf {gvcf_path}')

    j = b.new_job(f'ReblockGVCF', dict(sample=sample_name, project=project_name))
    j.image(resources.GATK_IMAGE)
    mem_gb = 8
    j.memory(f'{mem_gb}G')
    j.storage(f'50G')
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi',
        }
    )

    ref_fasta = resources.REF_FASTA
    ref_fai = resources.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )

    j.command(f"""
    export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
    gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

    function fail {{
      echo $1 >&2
      exit 1
    }}

    function retry {{
      local n=1
      local max=10
      local delay=30
      while true; do
        "$@" && break || {{
          if [[ $n -lt $max ]]; then
            ((n++))
            echo "Command failed. Attempt $n/$max:"
            sleep $delay;
          else
            fail "The command has failed after $n attempts."
          fi
        }}
      done
    }}

    GVCF=/io/batch/{sample_name}.g.vcf.gz
    GVCF_NODP=/io/batch/{sample_name}-nodp.g.vcf.gz
    REBLOCKED=/io/batch/{sample_name}-reblocked.g.vcf.gz

    # Retrying copying to avoid google bandwidth limits
    retry gsutil cp {gvcf_path} $GVCF
    retry gsutil cp {resources.NOALT_REGIONS} /io/batch/noalt-regions.bed

    # Copying reference data as well to avoid crazy logging costs 
    # for region requests
    retry gsutil cp {ref_fasta} /io/batch/{basename(ref_fasta)}
    retry gsutil cp {ref_fai}   /io/batch/{basename(ref_fai)}
    retry gsutil cp {ref_dict}  /io/batch/{basename(ref_dict)}

    # Reindexing just to make sure the index is not corrupted
    bcftools index --tbi /io/batch/{sample_name}.g.vcf.gz

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

    gatk --java-options "-Xms{mem_gb - 1}g" \\
    ReblockGVCF \\
    --reference /io/batch/{basename(ref_fasta)} \\
    -V $GVCF_NODP \\
    -do-qual-approx \\
    -O $REBLOCKED \\
    --create-output-variant-index true

    EXISTING_SN=$(bcftools view $REBLOCKED | awk '/^#CHROM/ {{ print $NF; exit }}')    
    bcftools view $REBLOCKED -T /io/batch/noalt-regions.bed \\
    | bcftools annotate -x INFO/DS \\
    | bcftools reheader -s <(echo "$EXISTING_SN {sample_name}") \\
    | bcftools view -Oz -o {j.output_gvcf['g.vcf.gz']}

    tabix -p vcf {j.output_gvcf['g.vcf.gz']}
    """)
    if output_path:
        b.write_output(j.output_gvcf, output_path.replace('.g.vcf.gz', ''))
    if depends_on:
        j.depends_on(*depends_on)
    return j
