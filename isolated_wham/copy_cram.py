#!/usr/bin/env python3

"""
runs a batch to acquire and store a CRAM file (+ index)
"""

from cpg_utils.config import config_retrieve, output_path
from cpg_utils.hail_batch import get_batch


job = get_batch().new_job("Copy Cram")
job.image(config_retrieve(['workflow', 'driver_image']))
job.storage('20Gi')

job.command(f'wget -O {job.cram} https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/1000Genomes/CRAM/HG00096/HG00096.cram')
job.command(f'wget -O {job.crai} https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/1000Genomes/CRAM/HG00096/HG00096.cram.crai')
get_batch().write_output(job.cram, output_path('HG00096.cram'))
get_batch().write_output(job.crai, output_path('HG00096.cram.crai'))
get_batch().run(wait=False)
