"""
Use Hail Batch to transfer VEP reference data cpg-reference bucket.
"""

import hailtop.batch as hb

from cpg_pipes import images, to_path
from cpg_pipes.hb.batch import setup_batch
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.providers.cpg import CpgStorageProvider
from cpg_pipes.refdata import RefData


MAKE_VEP_CACHE_TAR = False
MAKE_LOFTEE_TAR = False
PREPARE_MOUNTABLE_BUCKET = False


def main():
    """Entry point"""
    b = setup_batch('Copy VEP data')
    refs = RefData(CpgStorageProvider().get_ref_bucket())
    if MAKE_VEP_CACHE_TAR:
        _vep_cache(b, refs)
    if MAKE_LOFTEE_TAR:
        _loftee(b, refs)
    if PREPARE_MOUNTABLE_BUCKET:
        _uncompress(b, refs)
    _test(b, refs)
    res = b.run()
    res_status = res.status()
    assert res_status['state'] == 'success', str((res_status, res.debug_info()))


def _test(b: hb.Batch, refs: RefData):
    j = b.new_job('Test VEP mount')
    j.image(images.VEP_IMAGE)
    # gcsfuse works only with the root bucket, without prefix:
    base_bucket_name = refs.vep_bucket.drive
    data_mount = to_path(f'/{base_bucket_name}')
    j.cloudfuse(base_bucket_name, str(data_mount), read_only=True)
    vep_dir = data_mount / 'vep' / 'GRCh38'
    cmd = f"""
    ls {vep_dir}
    ls {vep_dir}/vep
    cat {vep_dir}/vep/homo_sapiens/105_GRCh38/info.txt
    """
    j.command(wrap_command(cmd))


def _uncompress(b: hb.Batch, refs: RefData):
    """
    Assuming tars are made and put on buckets, uncompresses them into a bucket
    to mount with gcsfuse.
    """
    j = b.new_job('Copy VEP cache with vep_install')
    j.image(images.VEP_IMAGE)
    j.storage(f'100G')

    cmd = f"""\
    mkdir /io/batch/uncompressed
    cd /io/batch/uncompressed
    df -h
    du -sh .

    tar -xvf {b.read_input(str(refs.vep_cache))}
    ls
    df -h
    du -sh .
    
    tar -xvf {b.read_input(str(refs.vep_loftee))}
    ls
    df -h
    du -sh .

    gsutil cp -r * {str(refs.vep_bucket)}/
    """
    j.command(wrap_command(cmd, setup_gcp=True))
    return j


def _vep_cache(b: hb.Batch, refs: RefData):
    """
    Prepare a tarball with VEP cache.
    """
    j = b.new_job('Copy VEP cache with vep_install')
    j.image(images.VEP_IMAGE)
    j.storage(f'30G')
    j.cpu(16)

    cmd = f"""\
    cd /io/batch
    
    # Because we are using biocontainers VEP image which, is conda-powered, 
    # we can use this conda-provided script:
    # https://bioconda.github.io/recipes/ensembl-vep/README.html#notes
    vep_install -a cf -s homo_sapiens -y GRCh38 -c vep --CONVERT
    tar -cvf {j.tar} vep
    """
    j.command(wrap_command(cmd))
    b.write_output(j.tar, refs.vep_cache)
    return j


def _loftee(b: hb.Batch, refs: RefData):
    """
    Prepare a tarball with LOFTEE ref data.
    """
    j = b.new_job('Prepare loftee reference bundle')
    j.image(images.VEP_IMAGE)
    j.storage(f'30G')
    cmd = f"""\
    cd /io/batch
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
    
    # loftee bigwig file (specific for grch38 loftee branch):
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

    tar -cvf {j.tar} .
    """
    j.command(wrap_command(cmd))
    b.write_output(j.tar, refs.vep_loftee)


if __name__ == '__main__':
    main()
