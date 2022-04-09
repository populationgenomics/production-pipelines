#!/usr/bin/env python3

"""
Jobs to find 65k somalier sites from gnomAD v3 using 
https://github.com/brentp/somalier/commit/a77fe13daeb4488d63c2819bde250fdd37ceb1cd
Problems with the currentl publicly shared somalier sites VCF:
 * is based of gnomAD 2.1.1 exomes, 
 * lifted over from GRCh37, 
 * and contains only 15k sites.
"""

import logging
from cpg_pipes import images, Namespace
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.providers.cpg import CpgStorageProvider

logger = logging.getLogger(__file__)

RESULT_VCF = CpgStorageProvider().get_ref_bucket() / 'somalier/v0/sites.hg38.vcf.gz'


pipe = create_pipeline(
    analysis_dataset='fewgenomes',
    name='find-somalier-sites',
    description='find 65k somalier sites',
    namespace=Namespace.MAIN,
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
gnomad_vcfs = [pipe.b.read_input(path) for path in gnomad_vcf_paths]
concat_j.command(
    wrap_command(
        f"""\
bcftools concat {" ".join(gnomad_vcfs)} -Oz -o {concat_j.gnomad_vcf}
""",
        monitor_space=True,
    )
)

make_sites_j = pipe.b.new_job('Make somalier sites')
make_sites_j.image(images.SOMALIER_IMAGE)
make_sites_j.storage('6T')
make_sites_j.cpu(4)
make_sites_j.command(
    wrap_command(
        f"""
cd /io/batch
somalier find-sites {concat_j.gnomad_vcf}
mv sites.vcf.gz {make_sites_j.sites_vcf}
"""
    )
)
pipe.b.write_output(make_sites_j.sites_vcf, RESULT_VCF)
