"""
Copy seqr reference data to gs://cpg-reference

Would first clean up existing data.
"""
import os

from cpg_pipes import to_path
from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.hb.command import wrap_command

src_seqr_bucket = 'gs://cpg-seqr-reference-data'
src_referen_path = f'{src_seqr_bucket}/GRCh38/all_reference_data/v2/combined_reference_data_grch38-2.0.3.ht'
src_clinvar_path = f'{src_seqr_bucket}/GRCh38/clinvar/clinvar.GRCh38.ht'

dst_seqr_bucket = 'gs://cpg-reference/seqr'


def main():
    b = setup_batch(
        description='Copy seqr reference data',
        billing_project=os.environ['HAIL_BILLING_PROJECT'],
        hail_bucket=to_path('gs://cpg-reference/hail-tmp'),
    )
    
    j = b.new_job('Copy reference data')
    j.image(os.environ['CPG_DRIVER_IMAGE'])
    cmd = f"""
    # gsutil -q rm -rf {dst_seqr_bucket}
    gsutil -q cp -r {src_referen_path} {dst_seqr_bucket}/
    gsutil -q cp -r {src_clinvar_path} {dst_seqr_bucket}/
    """
    j.command(wrap_command(cmd, setup_gcp=True))

    b.run(wait=False)


main()
