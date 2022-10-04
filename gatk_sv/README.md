# GATK-SV

Hail Batch driver for running the Broad's [GATK-SV workflow](https://github.com/broadinstitute/gatk-sv),
and importing results to Seqr, using the [cpg_utils/workflows library](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/workflows/README.md).

## Usage

To launch the workflow, run the `gatk_sv.py` script, and pass the following TOML configs directly as command line arguments. Replace `acute-care-test.toml` with another config if you want to process a different dataset.

```bash
python gatk_sv.py \
    ../configs/seqr.toml \
    ../configs/genome.toml \
    ../configs/acute-care-test.toml \
    ../configs/gatk-sv.toml
```

To run with the analysis runner, the usage is similar, but the configs are passed via the analysis runner command line arguments:

```bash
analysis-runner \
  --dataset acute-care \
  --description "GATK-SV" \
  --output-dir "gatk-sv" \
  --access-level test \
  --config ../configs/seqr.toml \
  --config ../configs/genome.toml \
  --config ../configs/acute-care-test.toml \
  --config ../configs/gatk-sv.toml \
  gatk_sv.py
```

The GATK-SV WDL workflow is run by a managed Cromwell instance via the analysis runner helper functions. The WDL workflow repository is [forked](https://github.com/populationgenomics/gatk-sv), and the latest commit to use is specified in `gatk_sv/gatk_sv.py` as a `GATK_SV_COMMIT` constant.

## Building reference and dockers

Running the workflow requires the Broad reference data copied to a `workflow/reference_prefix` bucket, and images copied to `workflow/image_registry_prefix`. In order to do that, use the following scripts:

```bash
bash sync_references.sh
python sync_images.py
```

The latter script will print a section for the `gatk-sv.toml` `[images]` section, that you'd need to insert there to make sure the workflow points to the latest tags. Also, when updating the GATK-SV fork, double check that the `[references.broad.sv.*]` sections in the config are up-to-date as well, and add new missing reference files if there are any.
