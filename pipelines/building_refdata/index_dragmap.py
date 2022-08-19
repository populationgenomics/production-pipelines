#!/usr/bin/env python3

"""
Create DRAGMAP reference index
"""

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes.providers.cpg.images import CpgImages
from cpg_pipes.providers.images import Images
from cpg_pipes.batch import setup_batch
from cpg_pipes.command import wrap_command
from cpg_pipes.providers.cpg.refdata import CpgRefData
from cpg_pipes.providers.refdata import RefData


def main():
    """
    Create index for DRAGMAP.
    """
    b = setup_batch('Create DRAGMAP index')
    refs = CpgRefData()
    images = CpgImages()
    j1 = _index_dragmap_job(b, refs, images)
    j2 = _test_dragmap_job(b, refs, images)
    j2.depends_on(j1)
    b.run(wait=False)


def _index_dragmap_job(b: hb.Batch, refs: RefData, images: Images) -> Job:
    """
    Creates the index for DRAGMAP
    """
    reference = refs.fasta_res_group(b)

    j = b.new_job('Index DRAGMAP')
    j.image(images.get('dragmap'))
    j.memory('standard')
    j.cpu(32)
    j.storage('40G')
    cmd = f"""\
    DIR=$(dirname {j.hash_table_cfg})

    dragen-os \\
    --build-hash-table true \\
    --ht-reference {reference.base} \\
    --ht-num-threads 32 \\
    --output-directory $DIR
    """
    j.command(wrap_command(cmd))
    for f in refs.dragmap_index_files:
        cmd += f'ln $DIR/{f} {getattr(j, f.replace(".", "_"))}\n'
    cmd += 'df -h; pwd; ls | grep -v proc | xargs du -sh'
    j.command(cmd)
    for f in refs.dragmap_index_files:
        b.write_output(getattr(j, f.replace('.', '_')), refs.dragmap_ref_bucket / f)
    return j


def _test_dragmap_job(b: hb.Batch, refs: RefData, images: Images) -> Job:
    dragmap_index = b.read_input_group(
        **{
            k.replace('.', '_'): refs.dragmap_ref_bucket / k
            for k in refs.dragmap_index_files
        }
    )
    fq1 = b.read_input('gs://cpg-seqr-test/batches/test/tmp_fq')
    j = b.new_job('Test DRAGMAP')
    j.image(images.get('dragmap'))
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
        wrap_command(
            f"""\
    dragen-os -r {dragmap_index} -1 {fq1} -2 {fq1} --RGID {sn} --RGSM {sn} |
    samtools sort -T $(dirname {j.sorted_bam})/samtools-sort-tmp -Obam -o {j.sorted_bam}
    """
        )
    )
    return j


if __name__ == '__main__':
    main()  # pylint: disable=E1120
