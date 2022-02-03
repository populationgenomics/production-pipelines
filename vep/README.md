# Preparing VEP for Hail query

Hail Query comes with VEP functionality, however it uses the image `konradjk/vep95_loftee:0.2` which is 3 years old and comes with an older VEP 95, while the lastest VEP to the date (27 Jan 2022) is 105. We want the new one because it annotates MANE transcripts. 

So firstly, we use a newer VEP image from biocontainers. We copy and modify [vep-GRCh38.sh](https://github.com/hail-is/hail/blob/main/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh) from the Hail codebase to set `VEP_DOCKER_IMAGE=quay.io/biocontainers/ensembl-vep:105.0--pl5262h4a94de4_0`.

Next, we need reference data cache matching VEP=105.

VEP cache bundle is huge, so we use Hail Batch to copy it from an FTP to the reference GCP bucket:

```bash
python copy_cache.py
```

For Loftee data, we can reuse Hail-provided bundle, which is already on a GCP bucket, so we can just copy it directly:

```bash
gsutil -u vlad-dev cp gs://hail-aus-sydney-vep/loftee-beta/GRCh38.tar \\
gs://cpg-reference/vep/loftee-beta.tar
```

Finally, copy the configs to the reference bucket:

```sh
gsutil cp vep105-GRCh38-loftee-gcloud.json vep-GRCh38.sh gs://cpg-reference/vep/
```
