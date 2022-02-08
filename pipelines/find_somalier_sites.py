#!/usr/bin/env python3

"""
Jobs to find 65k somalier sites from gnomAD v3 using 
https://github.com/brentp/somalier/commit/a77fe13daeb4488d63c2819bde250fdd37ceb1cd
Problems with the currentl publicly shared somalier sites VCF:
 * is based of gnomAD 2.1.1 exomes, 
 * lifted over from GRCh37, 
 * and contains only 15k sites.
"""

from os.path import join
import logging
from cpg_pipes import ref_data, images
from cpg_pipes.hailbatch import wrap_command
from cpg_pipes.pipeline import Pipeline

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)

RESULT_VCF = join(ref_data.REF_BUCKET, 'somalier/v0/sites.hg38.vcf.gz')


pipe = Pipeline(
    analysis_project='fewgenomes',
    name='find-somalier-sites',
    output_version='v0',
    namespace='main',
    title='find 65k somalier sites',
    check_smdb_seq_existence=False,
    keep_scratch=True,
)

concat_j = pipe.b.new_job('Make somalier sites')
concat_j.image(images.BCFTOOLS_IMAGE)
concat_j.storage('1000G')
concat_j.cpu(4)
gnomad_vcf_paths = [
    f'gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/'
    f'gnomad.genomes.v3.1.2.hgdp_tgp.chr{c}.vcf.bgz'
    for c in [str(i + 1) for i in range(22)] + ['X', 'Y']
]
gnomad_vcfs = [
    pipe.b.read_input(path) for path in gnomad_vcf_paths
]
concat_j.command(wrap_command(f"""\
bcftools concat {" ".join(gnomad_vcfs)} -Oz -o {concat_j.gnomad_vcf}
""", monitor_space=True))

make_sites_j = pipe.b.new_job('Make somalier sites')
make_sites_j.image(images.SOMALIER_IMAGE)
make_sites_j.storage('6T')
make_sites_j.cpu(4)
make_sites_j.command(wrap_command(f"""
cd /io/batch
somalier find-sites {concat_j.gnomad_vcf}
mv sites.vcf.gz {make_sites_j.sites_vcf}
"""))
pipe.b.write_output(make_sites_j.sites_vcf, RESULT_VCF)
