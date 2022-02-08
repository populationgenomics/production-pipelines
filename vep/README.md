# Preparing VEP for Hail Query

Hail Query comes with a `hl.vep()` function. However, it uses the image `konradjk/vep95_loftee:0.2`, which is built in 2019, and comes with an old VEP 95. While the lastest VEP to the date (27 Jan 2022) is 105. We want the latest VEP to date, because it annotates MANE transcripts, among just having more up-to-date data.

So firstly, we build an image with a newer VEP installed. See [../dockers/vep](../dockers/vep) for details.

Then we copy the dataproc [initialisation script](https://github.com/hail-is/hail/blob/cc0a051740f4de08408e6a2094ffcb1c3158ee9c/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh) from the Hail codebase and modify it to reflect the newly built image: [vep-GRCh38.sh](vep-GRCh38.sh)

Next, we rebuild the VEP cache, and LOFTEE reference data bundle to match VEP=v105. VEP cache bundle is huge, so we use Hail Batch to copy it from Ensemble FTP servers and the Broad Institute servers to the CPG reference GCP bucket.

```bash
python copy_cache.py
```

Then we copy the [VEP config](https://github.com/hail-is/hail/blob/cc0a051740f4de08408e6a2094ffcb1c3158ee9c/hail/python/hailtop/hailctl/hdinsight/resources/vep-GRCh38.json) from the Hail codebase, modify it to reflect the update reference bundle, and add `mane_select:String,mane_plus_clinical:String` into the schema:  [vep105-GRCh38-loftee-gcloud.json](vep105-GRCh38-loftee-gcloud.json).

Finally, we copy the initialization script and the config to the CPG reference bucket.

```sh
gsutil cp vep-GRCh38.sh vep105-GRCh38-loftee-gcloud.json gs://cpg-reference/vep/
```

Then you can start a cluster passing the initialisation script explicitly with `init`, instead of using the `vep` parameter:

```python
from analysis_runner import dataproc
j = dataproc.hail_dataproc_job(
    b,
    f'hail-query-script.py',
    max_age='4h',
    job_name='Run VEP',
    init=['gs://cpg-reference/vep/vep-GRCh38.sh'],
    worker_machine_type='n1-highmem-8',
    worker_boot_disk_size=200,
    secondary_worker_boot_disk_size=200,
    num_secondary_workers=20,
    num_workers=2,
)
```

You also need to make sure you set larger resources for the workers (the highmem machine type and larger storage). Hail would do that automatically with `--vep`, but with a custom --init we have to do that explicitly.

On that cluster, you can call the `hl.vep()` function with explicitly set `config` (otherwise it would attempt to get the VEP_CONFIG environment variable, which is only set with `--vep`):

```python
import hail as hl
mt = hl.vep(
    mt, 
    # We are not starting the cluster with --vep, instead passing custom
    # startup script with --init gs://cpg-reference/vep/vep-GRCh38.sh,
    # so VEP_CONFIG_URI will not be set, thus need to provide config
    # as a function parameter here:
    config='file:///vep_data/vep-gcloud.json'
)
```

Note that VEP would print errors messages that it can't find `phylocsf_data` and `gerp_exons` tables in the loftee DB:

```sh
DBD::SQLite::db prepare failed: no such table: phylocsf_data at /root/micromamba/share/ensembl-vep-105.0-0/LoF.pm line 565, <$fh> line 158458.
DBD::SQLite::db prepare failed: no such table: gerp_exons at /root/micromamba/share/ensembl-vep-105.0-0/gerp_dist.pl line 129, <$fh> line 158458.
```

It's expected as those tables are provided by `phylocsf_gerp.sql`, which are only provided for [GRCh37](https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/), but not for [GRCh38](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/). See [this issue](https://github.com/konradjk/loftee/issues/39).
