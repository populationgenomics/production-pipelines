"""
Run hap.py validation/concordance stats on a sample GVCF or joint VCF.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, image_path, fasta_res_group
from hailtop.batch import ResourceGroup

from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.targets import Sample


def happy(
    b,
    sample: Sample,
    vcf_or_gvcf: ResourceGroup,
    is_gvcf: bool,
    job_attrs: dict,
    output_path: Path | None = None,
):
    """
    Run hap.py validation/concordance stats on a sample GVCF or joint VCF.
    """
    if 'validation' not in get_config():
        logging.warning(
            'workflow/validation section is not defined in config, skip running hap.py'
        )
        return

    truth_sample_id = (
        get_config().get('validation', {}).get('sample_map').get(sample.participant_id)
    )
    if not truth_sample_id:
        return

    truth_vcf_path = reference_path(f'{truth_sample_id}/truth_vcf')
    truth_bed_path = reference_path(f'{truth_sample_id}/regions_bed')

    input_file: str
    if is_gvcf:
        input_file = vcf_or_gvcf['g.vcf.gz']
        prepy_params = '--convert-gvcf-to-vcf --filter-nonref'
        extract_sample_cmd = ''
    else:
        input_file = f'$BATCH_TMPDIR/{sample.id}.vcf.gz'
        prepy_params = ''
        # For multi-sample joint-called VCF, we need to extract our target sample.
        # we also want to strip FORMAT/AD and potentially other problematic fields,
        # which is done for gvcfs with `convert_gvcf_to_vcf` automatically:
        # https://github.com/Illumina/hap.py/blob/master/src/python/Tools/bcftools.py#L189
        # ...but not done for vcf.
        extract_sample_cmd = f"""
        bcftools view -s {sample.id} \
        {vcf_or_gvcf['vcf.gz']} -Ou | \
        bcftools annotate -x INFO,^FORMAT/GT,FORMAT/DP,FORMAT/GQ \
        -Oz -o {input_file}
        bcftools index {input_file}
        """

    # Calling regions
    seq_type = get_config()['workflow']['sequencing_type']
    eval_intervals_path = reference_path(f'broad/{seq_type}_evaluation_interval_lists')
    if seq_type == 'genome':
        # sparse regions, bcftools would loop through them
        regions_opt = '--restrict-regions'
    else:
        # dense regions, bcftools would use tabix to access each region in vcf
        regions_opt = '--target-regions'

    job_name = f'hap.py ({"GVCF" if is_gvcf else "VCF"})'
    job_attrs = (job_attrs or {}) | dict(tool='hap.py')
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('hap-py'))
    reference = fasta_res_group(b)
    res = STANDARD.set_resources(j, fraction=1)
    cmd = f"""\
    {extract_sample_cmd}
    
    grep -v ^@ {b.read_input(str(eval_intervals_path))} > intervals.bed
    head intervals.bed
    
    /opt/hap.py/bin/pre.py \
    --threads {res.get_nthreads()} \
    --pass-only \
    {prepy_params} \
    {input_file} \
    $BATCH_TMPDIR/pre-processed.vcf.gz \
    --reference {reference["base"]}

    # "--false-positives" means confident regions
    /opt/hap.py/bin/hap.py \
    --threads {res.get_nthreads()} \
    {b.read_input(str(truth_vcf_path))} \
    $BATCH_TMPDIR/pre-processed.vcf.gz \
    {regions_opt} intervals.bed \
    --false-positives {b.read_input(str(truth_bed_path))} \
    --report-prefix $BATCH_TMPDIR/prefix \
    --reference {reference["base"]}
    
    cp $BATCH_TMPDIR/prefix.summary.csv {j.summary_csv}
    """
    j.command(command(cmd))
    if output_path:
        b.write_output(j.summary_csv, str(output_path))
