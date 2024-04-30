# Run VEP in parallel using batch backend

This script runs VEP annotation on a VCF file using the Batch backend service. To run, a VCF file along with its tabix (.tbi) file are needed. Use conda to install the analysis-runner, then execute the following command:

```sh
analysis-runner \
    --config <toml config file> \
    --dataset <dataset> \
    --description <description>  \
    --output-dir <directory-to-output-annotated-ht> \
    --access-level <level> \
    vep_batch_backend.py \
    --vcf-path <path-to-vcf> \
    --output-ht <output-of-annotated-vcf>
```

As an example, the following invocation would run VEP annotation for a VCF file in the `test` bucket for the `tx-adapt` dataset.

```sh
analysis-runner \
    --config vep_config.toml \
    --dataset tx-adapt \
    --description "run VEP using batch backend"  \
    --output-dir "vep_batch/v0" \
    --access-level test \
    vep_batch_backend.py \
    # full path of VCF file
    --vcf-path "gs://cpg-tx-adapt-test/vep/indels.vcf.bgz" \
    --output-ht "vep/batch/annotated_indels.ht"
```
