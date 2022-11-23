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
from cpg_workflows.targets import Sample


def stripy(
    b,
    sample: Sample,
    cram_path: CramPath,
    target_loci: str,
    out_path: Path,
    log_path: Path,
    analysis_type: str = 'standard',
    write_to_bam: bool = False,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run STRipy

    Stripy is very inefficient at reading from cram files as it does not pass a local
    reference genome when opening the file reader. The right thing to do would be to patch this throughout
    the package and upstream. As a workaround we can either use the slow default (write_to_bam=False) or
    convert to a temporary bam then run stripy on that (write_to_bam=True) which is ~3x faster wall time but ~2x more
    expensive due to the local storage required.

    As implemented, this job will fail on a significant minority of our crams (aprox 1 in 5). No useful logs are produced
    even with the most verbose logging configuration set.
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

    reference = fasta_res_group(b)

    assert cram_path.index_path

    ## Atempt to mount bucket
    bucket = cram_path.path.drive
    print(f'bucket = {bucket}')
    bucket_mount_path = '/bucket'
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)

    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])
    mounted_cram_index_path = bucket_mount_path / '/'.join(cram_path.index_path.parts[2:])


    if write_to_bam:
        res = STANDARD.request_resources(storage_gb=100)
        res.set_to_job(j)
        write_bam_cmd = f"""
            BAM=$BATCH_TMPDIR/{cram_path.path.stem}.bam
            samtools view -b -@ 6 -T {reference.base}  $CRAM > $BAM
            samtools index $BAM
            ALIGNMENT=$BAM
        """
    else:
        res = STANDARD.request_resources(ncpu=2)
        res.set_to_job(j)
        write_bam_cmd = """ALIGNMENT=$CRAM"""

    if sample.pedigree.sex:
        sex_argument = f'--sex {str(sample.pedigree.sex).lower()}'
    else:
        sex_argument = ''

    cmd = f"""\

    ls -l {mounted_cram_path} {mounted_cram_index_path}

    cd ..
    git clone -b add-logging --single-branch https://gitlab.com/cassimons/stripy-pipeline.git stripy-test
    cd stripy-test
    chmod 755 batch.sh 
    
    # Increase logging to max verbosity. Needs to be passed as a config file so doing a quick an dirty edit
    # just edit of the default config on the fly and cat to the job log.
    sed 's/"log_flag_threshold": 1/"log_flag_threshold": -1/' config.json \
        > $BATCH_TMPDIR/config.json
    cat $BATCH_TMPDIR/config.json

    # CRAM=$BATCH_TMPDIR/{cram_path.path.name}
    # CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}
    CRAM={mounted_cram_path}
    CRAI={mounted_cram_index_path}

    # retry_gs_cp {str(cram_path.path)} $CRAM
    # retry_gs_cp {str(cram_path.index_path)} $CRAI

    {write_bam_cmd}
    
    python3 stri.py \\
        --genome hg38 \\
        --reference {reference.base} \\
        {sex_argument} \
        --output $BATCH_TMPDIR/ \\
        --input $ALIGNMENT \\
        --logflags {j.log_path} \\
        --config $BATCH_TMPDIR/config.json \\
        --analysis {analysis_type} \\
        --locus {target_loci}

    ls $BATCH_TMPDIR/
  
    cp $ALIGNMENT.html {j.out_path}

    """

    j.command(command(cmd, define_retry_function=True, monitor_space=True))
    b.write_output(j.log_path, str(log_path))
    b.write_output(j.out_path, str(out_path))

    return j
