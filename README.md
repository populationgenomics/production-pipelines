# Production pipelines

Hail Batch analysis pipelines used at CPG are located in the `pipelines` folder.

Example is the seqr loading pipeline, which takes CRAM/FASTQ data for a set of datasets (e.g. `acute-care` and `perth-neuro`), calls variants, annotates them, and creates an ElasticSearch index to use with seqr:

```sh
python pipelines/seqr_loader.py \
-n main \
--analysis-project seqr \
--input-project acute-care \
--input-project perth-neuro \
--output-project acute-care \
--ped-file gs://cpg-acute-care-main-upload/cpg_acute_positives_20211003_213917/acute-care-topup-mfranklin.ped \
--ped-file gs://cpg-acute-care-main-upload/acute-care-sm.ped \
--ped-file gs://cpg-perth-neuro-main-upload/perth-neuro-sm.ped \
--validate-smdb-analyses \
--check-smdb-seq-existence \
--check-intermediate-existence
--skip-sample CPG11783 \
--skip-sample CPG13326 \
--output-version v4-1 \
--keep-scratch
```

Another useful pipeline is `pipelines/pedigree.py` to verify inferred samples relationship and sex against a provided PED file(s):

```sh
python pipelines/pedigree.py \
-n main \
--analysis-project seqr \
--input-project acute-care \
--input-project perth-neuro \
--ped-file gs://cpg-acute-care-main-upload/acute-care-sm.ped \
--ped-file gs://cpg-perth-neuro-main-upload/perth-neuro-sm.ped \
--skip-sample CPG11783 \
--skip-sample CPG13326 \
--keep-scratch
```

The `cpg_pipes` package defines many handy python functions that which can be imported with `import cpg_pipes`. `cpg_pipes/jobs` defines functions that create Hail Batch Jobs for aligment, fastqc, deduplication, variant calling, VQSR, etc. For usage examples, see `pipelines/benchmarks/benchmark_alignment.py`, as well as other scripts in that folder.
