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

from cpg_utils.config import get_config, update_dict
from cpg_utils.hail_batch import reference_path, image_path

from workflows.command import wrap_command
from workflows.pipeline.pipeline import Pipeline

logger = logging.getLogger(__file__)


update_dict(
    get_config()['workflow'],
    {
        'name': 'find-somalier-sites',
        'description': 'find 65k somalier sites',
        'dataset': 'fewgenomes',
        'access_level': 'full',
        'keep_scratch': True,
    },
)
pipe = Pipeline()

results_vcf = reference_path('somalier_sites')


concat_j = pipe.b.new_job('Make somalier sites')
concat_j.image(image_path('bcftools'))
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
make_sites_j.image(image_path('somalier'))
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
pipe.b.write_output(make_sites_j.sites_vcf, str(results_vcf))
