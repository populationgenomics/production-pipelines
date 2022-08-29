#!/usr/bin/env python3

"""
Create DRAGMAP reference index
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils.hail_batch import command, image_path, reference_path, fasta_res_group
from cpg_utils.workflows.batch import get_batch

from jobs.align import DRAGMAP_INDEX_FILES


def main():
    """
    Create index for DRAGMAP.
    """
    b = get_batch('Create DRAGMAP index')
    j1 = _index_dragmap_job(b)
    j2 = _test_dragmap_job(b)
    j2.depends_on(j1)
    b.run(wait=False)


def _index_dragmap_job(b: hb.Batch) -> Job:
    """
    Creates the index for DRAGMAP
    """
    j = b.new_job('Index DRAGMAP')
    j.image(image_path('dragmap'))
    j.memory('standard')
    j.cpu(32)
    j.storage('40G')
    fasta_path = (
        reference_path('broad/dragmap_prefix') / 'Homo_sapiens_assembly38_masked.fasta'
    )
    cmd = f"""\
    DIR=$(dirname {j.hash_table_cfg})

    dragen-os \\
    --build-hash-table true \\
    --ht-reference {fasta_path} \\
    --ht-num-threads 32 \\
    --output-directory $DIR
    """
    j.command(command(cmd))
    for f in DRAGMAP_INDEX_FILES:
        cmd += f'ln $DIR/{f} {getattr(j, f.replace(".", "_"))}\n'
    cmd += 'df -h; pwd; ls | grep -v proc | xargs du -sh'
    j.command(cmd)
    for f in DRAGMAP_INDEX_FILES:
        b.write_output(
            getattr(j, f.replace('.', '_')), reference_path('broad/dragmap_prefix') / f
        )
    return j


def _test_dragmap_job(b: hb.Batch) -> Job:
    dragmap_index = b.read_input_group(
        **{
            k.replace('.', '_'): str(reference_path('broad/dragmap_prefix') / k)
            for k in DRAGMAP_INDEX_FILES
        }
    )
    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test DRAGMAP')
    j.image(image_path('dragmap'))
    j.cpu(32)
    j.memory('standard')
    j.storage('300G')
    j.declare_resource_group(
        output_cram={
            'cram': '{root}.cram',
            'crai': '{root}.crai',
        }
    )
    sn = 'TEST'
    j.command(
        command(
            f"""\
    dragen-os -r {dragmap_index} -1 {fq1} -2 {fq1} --RGID {sn} --RGSM {sn} |
    samtools sort -T $(dirname {j.sorted_bam})/samtools-sort-tmp -Obam -o {j.sorted_bam}
    """
        )
    )
    return j


if __name__ == '__main__':
    main()  # pylint: disable=E1120
