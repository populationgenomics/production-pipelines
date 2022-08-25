# Production workflows

# Seqr loader

A [Hail Batch](https://hail.is/docs/batch/service.html) and [metamist](https://github.com/populationgenomics/sample-metadata) driven workflows that loads genomics data into [Seqr](https://seqr.broadinstitute.org/).

## Usage

Set up a configuration TOML file, e.g. `configs/acute-care-test.toml`, and specify the datasets you want as inputs, sequencing type (genome or exome), and final stages you want to trigger, along with optional parameters:

```toml
[workflow]
datasets = ['acute-care', 'perth-neuro', 'validation']
# skip_samples = ['CPG11783', 'CPG255232', 'CPG255240']
sequencing_type = 'genome'
stages = ['MtToEs']
```

Set this file to `CPG_CONFIG_PATH` and run the analysis runner:

```bash
CPG_CONFIG_PATH=configs/acute-care-test.toml
analysis-runner \
  --dataset acute-care \
  --access-level test \
  -o seqr-loader \
  --description "test seqr loader" \
  main.py
```

### Installation

Requires Python 3.10. Clone the repository, and run:

```bash
pip install .
```
