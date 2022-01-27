# Preparing VEP for Hail query

Hail Query comes with VEP functionality, however it uses the image konradjk/vep95_loftee:0.2 which is 3 years old and comes with VEP v95, while the last one to the date (27 Jan 2022) is v105. We need the new one because it annotates MANE transcripts. So we first prepare a new image based on biocontainers (conda-powered) `ensemble-vep` image:

```bash
skopeo copy \
docker://quay.io/biocontainers/ensembl-vep:105.0--pl5262h4a94de4_0 \
docker://australia-southeast1-docker.pkg.dev/cpg-common/images/vep:105
```

For one-off jobs it's okay to use [conda-provided script](https://bioconda.github.io/recipes/ensembl-vep/README.html#notes): `vep_install -a cf -s homo_sapiens -y GRCh38 -c /vep --CONVERT`, however we want to keep the image small, and keep the reference data on a bucket (Hail does the same, so we can reuse their scripts; plus it makes it feasable to pulling the image on a local laptop).

VEP cache bundle is huge, so we use Hail Batch to copy it from an FTP to the reference GCP bucket:

```bash
python copy_cache.py
```

For Loftee data, we can reuse Hail-provided bundle, which is already on a GCP bucket, so we can just copy it directly:

```bash
gsutil -u vlad-dev cp gs://hail-aus-sydney-vep/loftee-beta/GRCh38.tar \\
gs://cpg-reference/vep/loftee-beta.tar
```
