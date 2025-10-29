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

    path_prefix = f'large_cohort/browser_data_download/v1-1/{data_type}'
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=path_prefix)

    paths = [f"gs://{bucket_name}/{b.name}" for b in blobs if b.name.endswith(".vcf.bgz")]

    b = get_batch()
    for path in paths:
        base = path.rsplit('/', 1)[1]
        output_path = path.replace('.vcf.gz', '.vcf.gz.csi')
        reindex_j = b.new_bash_job(f'Reindex browser VCF data download files for {path}')
        reindex_j.image(image_path('bcftools_121', '1.21-1'))
        reindex_j.storage('50Gi')
        reindex_j.command(
            f"""
set -euxo pipefail
IN="{base}"
gsutil -m cp "{path}" $BATCH_TMPDIR/
bcftools index -f -c "$BATCH_TMPDIR/$IN" -o ${reindex_j.out_idx}
""",
        )
        b.write_output(reindex_j.out_idx, output_path)

    b.run()


if __name__ == '__main__':
    init_batch()
    main()
