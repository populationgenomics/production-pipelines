#!/usr/bin/env python3

"""
Use Hail Batch to transfer VEP cache bundle into the cpg-reference bucket,
and prepare a GCP mount for cloudfuse.
"""

import click
import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import to_path
from cpg_utils.hail_batch import image_path, reference_path, command
from cpg_utils.workflows.batch import get_batch


@click.command()
@click.argument('vep_version')
def main(vep_version: str):
    """
    Build VEP VEP_VERSION cache and LOFTEE bundle
    """
    b = get_batch(f'Copy VEP v{vep_version} data')

    assert image_path('vep').split(':')[-1] == vep_version, (
        image_path('vep').split(':')[-1], vep_version
    )
    assert reference_path('vep_mount').parent.name == vep_version, (
        reference_path('vep_mount').parent.name == vep_version
    )

    j1 = _make_vep_cache_tar(b)
    j2 = _make_loftee_tar(b)
    j3 = _prepare_mountable_bucket(b)
    if j3:
        if j1:
            j3.depends_on(j1)
        if j2:
            j3.depends_on(j2)
    j4 = _test_mount(b)
    if j3:
        j4.depends_on(j3)

    res = b.run()
    res_status = res.status()
    assert res_status['state'] == 'success', str((res_status, res.debug_info()))


def _make_vep_cache_tar(b: hb.Batch) -> Job | None:
    """
    Prepare a tarball with VEP cache.
    """
    cache_path = reference_path('vep_mount').parent / 'cache.tar'
    if cache_path.exists():
        return None

    j = b.new_job('Copy VEP cache with vep_install')
    j.image(image_path('vep'))
    j.storage(f'30G')
    j.cpu(16)

    cmd = f"""\
    cd $BATCH_TMPDIR
    
    # Because we are using biocontainers VEP image which, is conda-powered, 
    # we can use this conda-provided script:
    # https://bioconda.github.io/recipes/ensembl-vep/README.html#notes
    # Use --NO_UPDATE to avoid user prompt if newer versions available
    vep_install -a cf -s homo_sapiens -y GRCh38 -c vep --CONVERT --NO_UPDATE
    tar -cvf {j.tar} vep
    """
    j.command(command(cmd))
    b.write_output(j.tar, str(cache_path))
    return j


def _make_loftee_tar(b: hb.Batch) -> Job | None:
    """
    Prepare a tarball with LOFTEE ref data.
    """
    loftee_path = reference_path('vep_mount').parent.parent / 'loftee.tar'
    if loftee_path.exists():
        return None

    j = b.new_job('Prepare LOFTEE bundle')
    j.image(image_path('vep'))
    j.storage(f'30G')
    cmd = f"""\
    cd $BATCH_TMPDIR
    mkdir loftee
    cd loftee
    
    # Following the LOFTEE docs https://github.com/konradjk/loftee/tree/grch38

    # conservation file:
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
    gunzip loftee.sql.gz
    
    # human_ancestor.fa file:
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
    
    # loftee bigwig file (specific for GRCh38 loftee branch):
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

    tar -cvf {j.tar} .
    """
    j.command(command(cmd))
    b.write_output(j.tar, str(loftee_path))
    return j


def _prepare_mountable_bucket(b: hb.Batch) -> Job | None:
    """
    Assuming tars are made and put on buckets, uncompresses them into a bucket
    to mount with gcsfuse.
    """
    if reference_path('vep_mount').exists():
        return None

    cache_path = reference_path('vep_mount').parent / 'cache.tar'
    loftee_path = reference_path('vep_mount').parent.parent / 'loftee.tar'

    j = b.new_job('Uncompresses bundles into a mount for cloudfuse')
    j.image(image_path('vep'))
    j.storage(f'100G')

    cmd = f"""\
    mkdir $BATCH_TMPDIR/uncompressed
    cd $BATCH_TMPDIR/uncompressed
    df -h
    du -sh .

    tar -xvf {b.read_input(str(cache_path))}
    ls
    df -h
    du -sh .
    
    tar -xvf {b.read_input(str(loftee_path))}
    ls
    df -h
    du -sh .

    gsutil cp -r * {str(reference_path('vep_mount'))}/
    """
    j.command(command(cmd, setup_gcp=True))
    return j


def _test_mount(b: hb.Batch) -> Job:
    j = b.new_job('Test VEP mount')
    j.image(image_path('vep'))
    # gcsfuse works only with the root bucket, without prefix:
    base_bucket_name = reference_path('vep_mount').drive
    data_mount = to_path(f'/{base_bucket_name}')
    j.cloudfuse(base_bucket_name, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(reference_path('vep_mount').parts[2:])
    cmd = f"""
    ls {vep_dir}
    ls {vep_dir}/vep
    cat {vep_dir}/vep/homo_sapiens/*/info.txt
    """
    j.command(command(cmd))
    return j


if __name__ == '__main__':
    main()
