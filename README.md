# Hail Batch Workflows

This repository provides implementations of genomics workflows using Hail Batch and the [cpg_utils/workflows library](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/workflows/README.md), specifically:

* WES/WGS Seqr Loader: FASTQ -> CRAM -> GVCF -> pVCF -> Matrix Table -> Elasticsearch index, with an ability to use parts of this pipeline as e.g. a single-sample germline variant calling workflow (FASTQ -> GVCF), joint-calling pipeline (GVCF -> pVCF), AS-VQSR pipeline, etc.
* (in progress) GATK-SV: CRAM -> VCF and more, orchestrates [GATK-SV WDL workflows](https://github.com/broadinstitute/gatk-sv) in Hail Batch.

## Seqr Loader

![uml](docs/seqr_loader.png)

Seqr Loader is a combination of the following [WARP](https://github.com/broadinstitute/warp)-inspired workflows:

* Single-sample whole genome and WES germline calling, FASTQ -> GVCF (BWA/DRAGMAP, GATK4 HaplotypeCaller)
* Re-alignment, CRAM -> GVCF (we use [BAZAM](https://github.com/ssadedin/bazam) to extract FASTQ from CRAM/BAM)
* Whole genome and WES joint-calling, GVCFs -> pVCF
* WES AS-VQSR and WGS VQSR workflows

As well as:

* [the Broad's seqr loading pipelines](https://github.com/broadinstitute/seqr-loading-pipelines)
* [Somalier](https://github.com/brentp/somalier) pedigree checks
* [MultiQC](https://github.com/ewels/MultiQC) QC reporting.

The workflows use Metamist as a source of FASTQ, CRAMs, and sample/participant metadata, and TOML configs for extended configuration (dataset and samples, Hail Batch parameters, Elasticsearch credentials, QC thresholds).

To run Seqr Loader, there is a script called `main.py`, which can take configs as input and should be submitted using the analysis-runner:

```bash
analysis-runner \
  --dataset seqr \
  --description "Seqr Loader" \
  --output-dir "seqr" \
  --access-level test \
  --config configs/seqr.toml \
  --config configs/genome.toml \
  --config configs/acute-care-test.toml
  main.py
```

Seqr Loader can be used partially, e.g. to just align and produce GVCFs (CRAM -> GVCF), do joint-calling (GVCF -> pVCF), run AS-VQSR, run matrix table annotation, create Elasticsearch index, create MultiQC reports, etc., if controlled by `workflows/first_stages`, `workflows/last_stages`, or `workflows/only_stages` parameters. It can also be followed by [the large cohort workflow](https://github.com/populationgenomics/large-cohort-pipeline) for:
	* VDS combiner to joint-call GVCFs directly into a VDS (sparse matrix table),
	* A dense matrix table subset for QC, PCA, kin, sex and ancestry inference,
	* PCA ancestry plots,
	* R markdown QC report.

## Stages and jobs

The `stages` folder provides implementations of the Seqr-Loader stages, that can be used altogether or in isolation. Stages use functions in the `jobs` module, that create jobs for different applications, e.g. alignment, deduplication, fastqc, haplotype caller, creating calling intervals, splitting VCF, running VQSR, etc. Those jobs can be called from the `queue_jobs()` method of a stage, or in isolation, as they don't need to know about the stage they are called from, but only need a `Batch` object. Feel free to explore both `stages` and `jobs` folders for examples and inspiration.
