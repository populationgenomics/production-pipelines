import click
from google.cloud import storage

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch, init_batch


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
    # list_blobs(prefix=...) returns EVERY object whose name starts with that prefix,
    # including those in the chrm/ subdirectory. We must exclude previously produced outputs.
    gvcf_files: list[str] = []
    for blob in blobs:
        name = blob.name
        # Skip anything in the output subdirectory
        if name.startswith(f'{prefix.rstrip("/")}/chrm/'):
            continue
        if not name.endswith('.gvcf.gz'):
            continue
        if 'chrM.hard' in name or 'chrM.hard-filtered' in name:
            continue
        gvcf_files.append(f'gs://{bucket_name}/{name}')

    gvcf_pairs = [(g, g + ".tbi") for g in gvcf_files]

    print(f"Found {len(gvcf_pairs)} GVCFs")

    for gvcf, tbi in gvcf_pairs:
        print(f"Processing {gvcf}")
        # create job
        j = b.new_bash_job(
            f'Subset {gvcf} tob-wgs gvcfs to chrM only',
            attributes={'tool': 'bcftools'},
        )
        j.image(image_path('bcftools'))
        j.storage('15Gi')

        # read in the gvcf files
        inputs = b.read_input_group(gvcf=gvcf, tbi=tbi)

        # declare output files
        j.declare_resource_group(
            output={
                'gvcf.gz': '{root}.chrM.hard-filtered.gvcf.gz',
                'gvcf.gz.csi': '{root}.chrM.hard-filtered.gvcf.gz.csi',
            },
        )

        # subset to chrM
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

        print(f'Writing to {output_chrm_gvcf_dir.rstrip("/")}/{gvcf.rsplit("/",1)[1].replace("hard","chrM.hard")}')

        b.write_output(
            j.output,
            f'{output_chrm_gvcf_dir.rstrip("/")}/{gvcf.rsplit("/",1)[1].replace("hard","chrM.hard")}',
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()
