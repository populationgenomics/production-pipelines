# Seqr loading pipeline

Seqr loading pipeline is essencially a WARP implementation on Hail Batch, followed by the Hail Query [seqr annotation pipeline](https://github.com/broadinstitute/hail-elasticsearch-pipelines):

* Aligning reads from a set of FASTQ pairs, or realigning from BAMs or CRAMs,
* Calling variants to produce per-sample GVCFs,
* Joint-calling,
* Loading data into a Hail MatrixTable,
* Annotating with VEP and other resources,
* Creating ElasticSearch index for each project.

To start the `seqr_loader` pipeline on 2 datasets `acute-care` and `perth-neuro` in the test namespace, using `seqr` as analysis dataset, and write the Elasticsearch index for `perth-neuro`:

```sh
python pipelines/seqr_loader.py \
-n main \
--analysis-dataset seqr \
--input-dataset acute-care \
--input-dataset perth-neuro \
--output-dataset perth-neuro \
--validate-smdb-analyses \
--check-smdb-seq-existence \
--check-intermediate-existence \
--skip-sample CPG11783 \
--skip-sample CPG13326 \
--output-version v1-0 \
--keep-scratch
```
