# Preparing VEP for Hail Query

Standard Hail Query `hl.vep()` function is limited to VEP version 95. Underneath it's using the image `konradjk/vep95_loftee:0.2` which is built in 2019 and hard pinned to VEP 95. 

To be able to run newer VEP, the following steps are required.

* Choose the VEP version you want to use. Make sure it's available on [Bioconda](https://anaconda.org/bioconda/ensembl-vep/files). Set the environment variable: `export VEP_VERSION=104.3`

* Add (if doing the first time) a Dockerfile into [images](https://github.com/populationgenomics/images/blob/07a2580c67886412ce1f0293274e7bd5e202a868/images/vep/Dockerfile), and trigger CI using the tag matching $VEP_VERSION.

* Add or update the `[images]` section in the [config template](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/config-template.toml#L146) to reflect this new version.

* Rebuild the VEP cache and LOFTEE reference data bundles. The VEP cache bundle is huge, so we use Hail Batch to copy it from Ensemble FTP servers and the Broad Institute servers to the CPG reference GCP bucket:

```bash
python copy-reference.py $VEP_VERSION
```

## Dataproc `hl.vep()`

Note that `hl.vep()` function works only with the spark backend on a Dataproc cluster, but is not supported on Batch backend.

* The dataproc initialisation script [dataproc-init.sh](dataproc-init.sh) is copied from the [Hail codebase](https://github.com/hail-is/hail/blob/cc0a051740f4de08408e6a2094ffcb1c3158ee9c/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh) and adjusted to pull the image from CPG artefact registry, and parameterised by VEP_VERSION.

Copy the new version of script into the bucket:

```shell
cat dataproc-init.sh | sed "s/__VEP_VERSION__/${VEP_VERSION}/g" | gsutil cp - gs://cpg-reference/vep/${VEP_VERSION}/dataproc/init.sh
```

* The JSON config script for dataproc [dataproc-config.json](dataproc-config.json) is also copied from the [Hail codebase](https://github.com/hail-is/hail/blob/cc0a051740f4de08408e6a2094ffcb1c3158ee9c/hail/python/hailtop/hailctl/hdinsight/resources/vep-GRCh38.json) and modified to reflect the newer VEP versions, and has `mane_select:String,mane_plus_clinical:String`.

The config is parametrised with `__VEP_VERSION__`, so pass it through `sed` when copying to the bucket:

```sh
cat dataproc-config.json | sed "s/__VEP_VERSION__/${VEP_VERSION}/g" | gsutil cp - gs://cpg-reference/vep/${VEP_VERSION}/dataproc/config.json
```

* After all set up, you can start a Dataproc cluster passing the initialization script explicitly with `init`, instead of using the `vep` parameter:

```python
from analysis_runner import dataproc
from cpg_utils.workflows.batch import get_batch
vep_version = ...
j = dataproc.hail_dataproc_job(
    get_batch('Run VEP'),
    f'some-hail-query-script.py',
    max_age='4h',
    job_name=f'Run VEP {vep_version}',
    init=[f'gs://cpg-reference/vep/GRCh38-{vep_version}-dataproc.sh'],
    worker_machine_type='n1-highmem-8',
    worker_boot_disk_size=200,
    secondary_worker_boot_disk_size=200,
    num_secondary_workers=20,
    num_workers=2,
)
```

Make sure you set larger resources for the workers (the highmem machine type and larger storage). Hail would do that automatically with `--vep`, but with a custom `--init` we have to do that explicitly.

On that cluster, you can call the `hl.vep()` function with explicitly set `config` (otherwise it would attempt to get the VEP_CONFIG environment variable, which is only set with `--vep`):

```python
import hail as hl
mt = ...
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
