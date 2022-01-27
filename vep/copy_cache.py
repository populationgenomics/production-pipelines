"""
Use Hail Batch to transfer VEP cache to the cpg-reference bucket
"""
from cpg_pipes import resources
from cpg_pipes.hailbatch import setup_batch, wrap_command

b = setup_batch('Copy VEP cache')
j = b.new_job('Copy VEP cache')
j.cpu(16)
j.image(resources.DRIVER_IMAGE)

cmd = """\
curl http://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz | gsutil cp - gs://cpg-reference/vep/homo_sapiens_vep_105_GRCh38.tar.gz\
"""

j.command(wrap_command(cmd, setup_gcp=True))
b.run(wait=False)
