import click
from google.cloud import storage

import hail as hl

from cpg_utils import Path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch, image_path, init_batch


@click.command()
@click.option('--bucket', type=str, required=True, help='GCS bucket name where input gvcfs are stored.')
@click.option('--path-prefix', type=str, required=True, help='Path prefix within the bucket to look for gvcfs.')
@click.option('--output-chrm-gvcf-dir', type=str, required=True, help='GCS directory to write chrM gvcfs to.')
def main(bucket: str, path_prefix: str, output_chrm_gvcf_dir: str):
    """
    Subset tob gvcfs to chrM only.
    :param input_gvcf_dir: Path to input gvcf directory.
    :param output_chrm_gvcf_dir: Path to output chrm gvcf directory.
    """
    init_batch()
    b = get_batch()

    client = storage.Client()
    bucket_name = bucket
    prefix = path_prefix

    blobs = client.list_blobs(bucket_name, prefix=prefix)

    gvcf_files = [f"gs://{bucket_name}/{b.name}" for b in blobs if b.name.endswith(".gvcf.gz")]

    gvcf_pairs = [(g, g + ".tbi") for g in gvcf_files]

    print(f"Found {len(gvcf_pairs)} GVCFs")

    for gvcf, tbi in gvcf_pairs:
        # create job
        j = b.new_bash_job(
            f'Subset {gvcf} tob-wgs gvcfs to chrM only',
            attributes={'tool': 'bcftools'},
        )
        j.image(image_path('bcftools'))

        j.declare_resource_group(
            output={
                'gvcf.gz': '{root}.chrM.hard.gvcf.gz',
                'gvcf.gz.tbi': '{root}.chrM.hard.gvcf.gz.tbi',
            },
        )

        # read in the gvcf files
        inputs = b.read_input_group(gvcf=gvcf, tbi=tbi)

        j.command(
            f"""
            set -euxo pipefail

            IN_VCF="{inputs['gvcf']}"

            OUT_VCF="{j.output['gvcf.gz']}"

            bcftools view -r chrM -Oz -o "$OUT_VCF" "$IN_VCF"
            bcftools index "$OUT_VCF"

            ls -l $(dirname "$OUT_VCF")
            """,
        )
        b.write_output(
            j.output['gvcf.gz'],
            f'{output_chrm_gvcf_dir.rstrip("/")}/{gvcf.rsplit("/",1)[1].replace("hard","chrM.hard")}',
        )
        b.write_output(
            j.output['gvcf.gz.tbi'],
            f'{output_chrm_gvcf_dir.rstrip("/")}/{gvcf.rsplit("/",1)[1].replace("hard","chrM.hard")}.tbi',
        )


if __name__ == '__main__':
    main()
