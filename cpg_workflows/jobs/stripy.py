"""
Create Hail Batch jobs to run STRipy
"""

from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse
from cpg_workflows.targets import Sample


def stripy(
    b,
    sample: Sample,
    cram_path: CramPath,
    target_loci: str,
    out_path: Path,
    log_path: Path,
    json_path: Path,
    analysis_type: str = 'standard',
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

    reference = fasta_res_group(b)

    # Stripy accesses a relatively small number of discrete regions from each cram
    # accessing the cram via cloudfuse is faster than localising the full cram   
    bucket = cram_path.path.drive
    print(f'bucket = {bucket}')
    bucket_mount_path = to_path('/bucket')
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)
    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    if sample.pedigree.sex:
        sex_argument = f'--sex {str(sample.pedigree.sex).lower()}'
    else:
        sex_argument = ''

    cmd = f"""\
    # Increase logging to max verbosity and output json results. Needs to be passed as a config file so doing a 
    # quick an dirty edit of the default config on the fly and cat to the job log.
    sed 's/"log_flag_threshold": 1/"log_flag_threshold": -1/' config.json \
        | sed 's/"output_json": false/"output_json": true/' \
        > $BATCH_TMPDIR/config.json
    cat $BATCH_TMPDIR/config.json

    python3 stri.py \\
        --genome hg38 \\
        --reference {reference.base} \\
        {sex_argument} \
        --output $BATCH_TMPDIR/ \\
        --input {mounted_cram_path} \\
        --logflags {j.log_path} \\
        --config $BATCH_TMPDIR/config.json \\
        --analysis {analysis_type} \\
        --locus {target_loci}

    ls $BATCH_TMPDIR/
  
    cp $BATCH_TMPDIR/{cram_path.path.name}.html {j.out_path}
    cp $BATCH_TMPDIR/{cram_path.path.name}.json {j.json_path}

    """

    j.command(command(cmd))
    b.write_output(j.log_path, str(log_path))
    b.write_output(j.json_path, str(json_path))
    b.write_output(j.out_path, str(out_path))

    return j
