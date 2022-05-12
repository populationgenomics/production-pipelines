"""
Copy seqr reference data to gs://cpg-reference

Would first clean up existing data.
"""
import os
from os.path import basename

from cpg_pipes import to_path
from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.hb.command import wrap_command, python_command, GCLOUD_CMD

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
    python_cmd = f"""
    import logging
    logger = logging.getLogger(__file__)
    logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
    logger.setLevel(logging.INFO)
    
    import asyncio
    import hail as hl
    asyncio.get_event_loop().run_until_complete(
        hl.init_batch(
            default_reference='GRCh38',
            billing_project='{os.environ["HAIL_BILLING_PROJECT"]}',
            remote_tmpdir='{to_path("gs://cpg-reference/hail-tmp")}',
        )
    )
    ht = hl.read_table('{src_referen_path}')
    ht.write('{dst_seqr_bucket}/{basename(src_referen_path)}', overwrite=True)

    ht = hl.read_table('{src_clinvar_path}')
    ht.write('{dst_seqr_bucket}/{basename(src_clinvar_path)}', overwrite=True)
    """

    cmd = f"""
    set -o pipefail
    set -ex
    {GCLOUD_CMD}

    cat << EOT >> script.py
    {python_cmd}
    EOT
    python3 script.py
    """
    j.command(wrap_command(cmd, monitor_space=True))
    b.run(wait=False)


main()
