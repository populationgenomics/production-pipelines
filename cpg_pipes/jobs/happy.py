"""
Run hap.py validation/concordance stats on a sample GVCF or joint VCF.
"""

import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, image_path, fasta_res_group
from hailtop.batch import ResourceGroup

from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.targets import Sample

logger = logging.getLogger(__file__)


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
    if 'validation' not in get_config()['workflow']:
        logger.warning(
            'workflow/validation section is not defined in config, skip running hap.py'
        )
        return

    sequencing_type = get_config()['workflow']['sequencing_type']
    if sequencing_type not in get_config()['workflow']['validation']:
        logger.warning(
            f'workflow/validation/{sequencing_type} section is not defined in config, skip running hap.py'
        )
        return

    truth_sample_id = get_config()['workflow']['validation'].get(sample.external_id)
    if not truth_sample_id:
        return

    truth_vcf_path = reference_path(
        f'validation/{sequencing_type}/{truth_sample_id}/truth'
    )
    truth_bed_path = reference_path(
        f'validation/{sequencing_type}/{truth_sample_id}/regions'
    )
    input_file: str

    if is_gvcf:
        input_file = vcf_or_gvcf['g.vcf.gz']
        happy_params = '--convert-gvcf-to-vcf --filter-nonref'
        extract_sample_cmd = ''
    else:
        input_file = f'/io/{sample.id}.vcf.gz'
        happy_params = ''
        # For multi-sample joint-called VCF, we need to extract our target sample
        extract_sample_cmd = f"""
        bcftools view -s {sample.id} \
        {vcf_or_gvcf['vcf.gz']} \
        -Oz -o {input_file}
        bcftools index {input_file}
        """

    jname = 'Happy'
    job_attrs = (job_attrs or {}) | dict(tool='hap.py')
    j = b.new_job(jname, job_attrs)
    j.image(image_path('happy'))
    reference = fasta_res_group(b)
    res = STANDARD.set_resources(j, fraction=1)
    cmd = f"""\
    {extract_sample_cmd}
    
    /opt/hap.py/bin/pre.py \
    --threads {res.get_nthreads()} \
    --pass-only \
    {happy_params} \
    {input_file} \
    /io/pre-processed.vcf.gz \
    -r {reference["base"]}
    
    /opt/hap.py/bin/hap.py \
    --threads {res.get_nthreads()} \
    {b.read_input(str(truth_vcf_path))} \
    /io/pre-processed.vcf.gz \
    -f {b.read_input(str(truth_bed_path))} \
    -o /io/batch/prefix \
    -r {reference["base"]}
    
    cp /io/batch/prefix.summary.csv {j.summary_csv}
    """
    j.command(wrap_command(cmd))
    if output_path:
        b.write_output(j.summary_csv, str(output_path))
