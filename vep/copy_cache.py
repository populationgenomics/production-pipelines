"""
Use Hail Batch to transfer VEP cache to the cpg-reference bucket.

Because we are using biocontainers VEP image which, is conda-powered, we can use 
this [conda-provided script](https://bioconda.github.io/recipes/ensembl-vep/README.html#notes):
`vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep --CONVERT` to pull VEP's cache.
"""

from cpg_pipes.hailbatch import setup_batch, wrap_command

b = setup_batch('Copy VEP cache with vep_install')
j = b.new_job('Copy VEP cache')
j.image('quay.io/biocontainers/ensembl-vep:105.0--pl5262h4a94de4_0')
j.storage(f'30G')
j.cpu(16)

cmd = f"""\
cd /io/batch
vep_install -a cf -s homo_sapiens -y GRCh38 -c vep --CONVERT
tar -cvf {j.tar} vep
"""

j.command(wrap_command(cmd))
b.write_output(j.tar, 'gs://cpg-reference/vep/homo_sapiens_vep_105_GRCh38.tar')
b.run(wait=False)
