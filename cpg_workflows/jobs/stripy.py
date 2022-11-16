"""
Create Hail Batch jobs to run STRipy
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import image_path, fasta_res_group, reference_path
from cpg_utils.hail_batch import command
from cpg_workflows.resources import HIGHMEM, STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse, exists


def stripy(
    b,
    cram_path: CramPath,
    out_pdf_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run STRipy
    """
    if can_reuse(
        [
            out_pdf_path,
        ],
        overwrite,
    ):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'stripy'}
    j = b.new_job('STRipy', job_attrs)
    j.image(image_path('stripy'))
    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)
    reference = fasta_res_group(b)

    assert cram_path.index_path
    cmd = f"""\
    CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

    # Retrying copying to avoid google bandwidth limits
    retry_gs_cp {str(cram_path.path)} $CRAM
    retry_gs_cp {str(cram_path.index_path)} $CRAI

    python3 stri.py \\
        --genome hg38 \\
        --reference {reference.base} \\
        --output $BATCH_TMPDIR/ \\
        --locus AFF2,ATXN3,HTT,PHOX2B \\
        --input $CRAM

    ls $BATCH_TMPDIR/
    cp $BATCH_TMPDIR/foo {j.out_pdf_path}
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.out_pdf_path, str(out_pdf_path))
    return j
