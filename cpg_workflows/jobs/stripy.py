"""
Create Hail Batch jobs to run STRipy
"""

import json

import hailtop.batch as hb

from cpg_utils import Path, to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import command, fasta_res_group
from cpg_workflows.filetypes import CramPath
from cpg_workflows.resources import STANDARD
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import can_reuse


def stripy(
    b: hb.batch.Batch,
    sequencing_group: SequencingGroup,
    cram_path: CramPath,
    target_loci: str,
    out_path: Path,
    log_path: Path,
    json_path: Path,
    custom_loci_path: str = '',
    analysis_type: str = 'standard',
    stripy_config: dict | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> hb.batch.job.Job | None:
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

    config_path = 'config.json'
    if stripy_config:
        j.command("echo original config:")
        j.command(f"cat {config_path}")
        j.command(
            f"echo $(cat {config_path} | jq '. * $p' {config_path} --argjson p '{json.dumps(stripy_config)}') > $BATCH_TMPDIR/config_updated.json",
        )
        config_path = '$BATCH_TMPDIR/config_updated.json'

    reference = fasta_res_group(b)

    # Stripy accesses a relatively small number of discrete regions from each cram
    # accessing the cram via cloudfuse is faster than localising the full cram
    bucket = cram_path.path.drive
    print(f'bucket = {bucket}')
    bucket_mount_path = to_path('/bucket')
    j.cloudfuse(bucket, str(bucket_mount_path), read_only=True)
    mounted_cram_path = bucket_mount_path / '/'.join(cram_path.path.parts[2:])
    assert cram_path.index_path  # keep mypy happy as index_path is optional
    mounted_cram_index_path = bucket_mount_path / '/'.join(cram_path.index_path.parts[2:])

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    if sequencing_group.pedigree.sex and str(sequencing_group.pedigree.sex).lower() != 'unknown':
        sex_argument = f'--sex {str(sequencing_group.pedigree.sex).lower()}'
    else:
        sex_argument = ''

    if custom_loci_path:
        custom_loci_input = b.read_input(str(custom_loci_path))
        custom_loci_argument = f'--custom {custom_loci_input}'
    else:
        custom_loci_argument = ''

    cmd = f"""\
    cat {config_path}

    ln -s {mounted_cram_path} {sequencing_group.id}__{sequencing_group.external_id}.cram
    ln -s {mounted_cram_index_path} {sequencing_group.id}__{sequencing_group.external_id}.crai

    python3 stri.py \\
        --genome hg38 \\
        --reference {reference.base} \\
        {sex_argument} \
        --output $BATCH_TMPDIR/ \\
        --input {sequencing_group.id}__{sequencing_group.external_id}.cram  \\
        --logflags {j.log_path} \\
        --config {config_path} \\
        --analysis {analysis_type} \\
        --locus {target_loci} \\
        {custom_loci_argument}


    ls $BATCH_TMPDIR/

    cp $BATCH_TMPDIR/{sequencing_group.id}__{sequencing_group.external_id}.cram.html {j.out_path}
    cp $BATCH_TMPDIR/{sequencing_group.id}__{sequencing_group.external_id}.cram.json {j.json_path}

    """

    j.command(command(cmd))
    b.write_output(j.log_path, str(log_path))
    b.write_output(j.json_path, str(json_path))
    b.write_output(j.out_path, str(out_path))

    return j
