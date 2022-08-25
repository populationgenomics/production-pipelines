# Production workflows

# Seqr loader

A [Hail Batch](https://hail.is/docs/batch/service.html) and [metamist](https://github.com/populationgenomics/sample-metadata) driven workflows that loads genomics data into [Seqr](https://seqr.broadinstitute.org/).

## Usage

Prepare a configuration TOML file, e.g. `configs/seqr-test.toml` or `configs/seqr.toml`, where you can specify the datasets you want as inputs, sequencing type (genome or exome), and configure target stages of the pipeline. You can use multiple configs.

```bash
analysis-runner \
  --dataset acute-care \
  --access-level test \
  --output-dir "seqr-loader" \
  --description "test seqr loader" \
  --config configs/seqr-test.toml \
  --config configs/genome.toml \
  main.py
```
