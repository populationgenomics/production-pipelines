import os

import click
from google.cloud import storage

import hail as hl

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch


@click.command()
@click.option(
    '--data-type',
    type=click.Choice(['exomes', 'genomes', 'joint'], case_sensitive=False),
    required=True,
    help='Dataset type.',
)
@click.option(
    '--bucket-name',
    type=str,
    help='GCS bucket name where browser VCF data download files are stored.',
)
def main(data_type: str, bucket_name: str):
    """
    Current VCF files are unable to be viewed correctly. This is due to the indexing files.
    This script reindexes browser VCF data download files.
    """

    PREFIX_ROOT = {
        'exomes': 'exome/large_cohort/browser_data_download/v1-1/exomes',
        'genomes': 'large_cohort/browser_data_download/v1-1/genomes',
        'joint': 'large_cohort/browser_data_download/v1-1/joint',
    }

    path_prefix = PREFIX_ROOT[data_type]
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=path_prefix)

    paths = [f"gs://{bucket_name}/{b.name}" for b in blobs if b.name.endswith(".vcf.bgz")]

    b = get_batch()
    for path in paths:
        input_file = b.read_input(path)
        output_path = path.replace('.vcf.bgz', '.vcf.bgz.csi')
        reindex_j = b.new_bash_job(f'Reindex browser VCF data download files for {path}')
        reindex_j.image(image_path('bcftools_121', '1.21-1'))
        reindex_j.storage('50Gi')
        reindex_j.command(
            f"""
set -euxo pipefail
IN="{input_file}"
bcftools index -f -c "$IN" -o {reindex_j.out_idx}
""",
        )
        b.write_output(reindex_j.out_idx, output_path)

    b.run()


if __name__ == '__main__':
    init_batch()
    main()
