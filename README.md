# Production pipelines

Helper functions to run Hail Batch pipelines at the CPG.

* `cpg_pipes` contains the codebase for the `cpg-pipes` pip package, which includes multiple utility functions:
  * `cpg_pipes/jobs`: create Job objects for various bioinformatics purposes,
  * `cpg_pipes/hb`: interact with Hail Batch and build new jobs,
  * `cpg_pipes/buckets`: interact with objects on buckets (GCP, and Azure in the future),
  * `cpg_pipes/images`: Docker image tags hosted by the CPG,
  * `cpg_pipes/ref_data`: bioinformatics reference data hosted by the CPG,
* `batches`: example Batches that use jobs from `cpg_pipes/jobs`

Batches must be run with the [analysis-runner](https://github.com/populationgenomics/analysis-runner):

```bash
# make sure you've installed https://anaconda.org/cpg/analysis-runner
analysis-runner \
  --access-level test \
  --dataset fewgenomes \
  --description "Run VEP" \
  --output-dir "$(whoami)-test-vep" \
  batches/vep.py \
  --vcf gs://cpg-fewgenomes-test/unittest/inputs/chr20/genotypegvcfs/vqsr.vcf.gz \
  --project fewgenomes
```
