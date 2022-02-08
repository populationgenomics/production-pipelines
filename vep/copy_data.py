"""
Use Hail Batch to transfer VEP reference data cpg-reference bucket.
"""

from cpg_pipes.hailbatch import setup_batch, wrap_command, resources


def _vep_cache(b):
    j = b.new_job('Copy VEP cache with vep_install')
    j.image(resources.VEP_IMAGE)
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
    b.write_output(j.tar, 'gs://cpg-reference/vep/homo_sapiens_vep_105_GRCh38.tar')
    return j


def _loftee(b):
    j = b.new_job('Prepare loftee reference bundle')
    j.image(resources.VEP_IMAGE)
    j.storage(f'30G')
    cmd = f"""\
    cd /io/batch
    mkdir loftee
    cd loftee
    
    # Following the LOFTEE docs https://github.com/konradjk/loftee/tree/grch38

    # conservation file:
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
    
    # human_ancestor.fa file:
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
    
    # loftee bigwig file (specific for grch38 loftee branch):
    wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

    tar -cvf {j.tar} .
    """
    j.command(wrap_command(cmd))
    b.write_output(j.tar, 'gs://cpg-reference/vep/loftee_GRCh38.tar')


b = setup_batch('Copy VEP data')
# _vep_cache(b)
_loftee(b)
b.run(wait=False)
