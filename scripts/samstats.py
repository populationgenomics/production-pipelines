#!/usr/bin/env python3

"""
Runs sam stats on one dumb file
"""

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch
from cpg_utils.hail_batch import command, fasta_res_group
from cpg_utils.config import image_path


job = get_batch().new_job("Copy Cram")
job.image(image_path('samtools'))

job.storage('20Gi')

reference = fasta_res_group(b)

cmd = f"""\
    CRAM=$BATCH_TMPDIR/in.cram
    CRAI=$BATCH_TMPDIR/in.cram.crai

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp gs://cpg-seqr-test/thousand_genomes_copy/HG00096.cram $CRAM
    retry_gs_cp gs://cpg-seqr-test/thousand_genomes_copy/HG00096.cram.crai $CRAI

    samtools stats \\
    --reference {reference.base} \\
    $CRAM > {job.output_stats}
    """

job.command(command(cmd, define_retry_function=True))
get_batch().write_output(job.output_stats, 'gs://cpg-seqr-test/thousand_genomes_copy/HG00096.cram.samstats')

get_batch().run(wait=False)

job.command(command(cmd, define_retry_function=True))
