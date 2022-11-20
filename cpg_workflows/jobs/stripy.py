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
    out_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run STRipy
    """
    if can_reuse(
        [
            out_path,
        ],
        overwrite,
    ):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'stripy'}
    j = b.new_job('STRipy', job_attrs)
    j.image(image_path('stripy'))
    # res = STANDARD.request_resources(ncpu=2, storage_gb=100)
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

    # samtools view -b -@ 2 -T {reference.base}  $CRAM > $CRAM.bam
    # samtools index $CRAM.bam

    # --analysis extended 
    python3 stri.py \\
        --genome hg38 \\
        --reference {reference.base} \\
        --output $BATCH_TMPDIR/ \\
        --input $CRAM \\
        --locus AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,RAPGEF2,RFC1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3

    echo "BATCH_TMPDIR = $BATCH_TMPDIR"
    ls $BATCH_TMPDIR/
  
    # cp $CRAM.bam.html {j.out_path}
    cp $CRAM.html {j.out_path}

    """

    j.command(command(cmd, define_retry_function=True, monitor_space=True))
    b.write_output(j.out_path, str(out_path))
    return j
