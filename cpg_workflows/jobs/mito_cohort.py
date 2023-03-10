"""
Create Hail Batch jobs to call mitochondrial SNVs
"""
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse

from cpg_workflows.jobs import picard
from cpg_workflows.mito_pipeline_scripts import annotate_coverage as annotate_coverage_script

def split_multi_allelics(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    remove_non_pass_sites: bool = False,
    job_attrs: dict | None = None,
) -> Job:
    """
    Splits multi allelics and removes non pass sites
    Uses LeftAlignAndTrimVariants to split then optionally use SelectVariants to select
    only passing variants.
    Args:
        vcf: Input vcf file.
        reference: chrM reference fasta.
        remove_non_pass_sites:
    Output:
        output_vcf: Final vcf file.
    Cmd from:
    https://github.com/broadinstitute/gatk/blob/4ba4ab5900d88da1fcf62615aa038e5806248780/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L600
    """
    job_attrs = job_attrs or {}
    j = b.new_job('split_multi_allelics', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        split_vcf={'vcf.gz': '{root}.vcf.gz'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz'}
    )

    cmd = f"""
        gatk LeftAlignAndTrimVariants \
            -R {reference.base} \
            -V {vcf['vcf.gz']} \
            -O {j.split_vcf['vcf.gz']} \
            --split-multi-allelics \
            --dont-trim-alleles \
            --keep-original-ac
        """
    if remove_non_pass_sites:
        cmd += f"""
            gatk SelectVariants \
                -V {j.split_vcf['vcf.gz']} \
                -O {j.output_vcf['vcf.gz']} \
                --exclude-filtered
        """
    else:
        cmd += f"""
            mv {j.split_vcf['vcf.gz']} {j.output_vcf['vcf.gz']}
        """

    j.command(command(cmd, define_retry_function=True))

    return j


def annotate_coverage(
    b,
    base_level_coverage_by_sid: dict[str, Path],
    out_mt: hb.ResourceFile,
    # haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Combine individual mitochondria coverage files and outputs a hail table with coverage annotations
    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_annotate_coverage', job_attrs)
    # j.image(image_path('haplocheckcli'))
    # anno_job.image(get_config()['workflows']['driver_image'])

    res = STANDARD.request_resources(ncpu=8)
    res.set_to_job(j)

    # generate tsv string to use as input file
    # required format: "tsv of participant_id, base_level_coverage_metrics, sample"
    tsv_string = ""
    for sid, path in base_level_coverage_by_sid.items():
        tsv_string += f'{sid}\t{path}\t{sid}\n'

    cmd = f"""
        # build input file:
        echo "{tsv_string}" > input.tsv

        # Run query job
        python {annotate_coverage_script.__file__} \
            --input-tsv input.tsv \
            --output-ht {j.out_mt}
            --temp-dir $BATCH_TMPDIR/mt
        """

    j.command(command(cmd, setup_gcp=True, monitor_space=True))
    b.write_output(j.out_mt, str(out_mt))

    return j
