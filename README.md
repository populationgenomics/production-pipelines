# Hail Batch Workflows

This repository provides implementations of genomics workflows using Hail Batch and the [cpg_utils/workflows library](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/workflows/README.md), specifically:

* WES/WGS Seqr Loader: FASTQ -> CRAM -> GVCF -> pVCF -> Matrix Table -> Elasticsearch index, with an ability to use parts of this pipeline as e.g. a single-sample germline variant calling workflow (FASTQ -> GVCF), joint-calling pipeline (GVCF -> pVCF), AS-VQSR pipeline, etc.
* (in progress) GATK-SV: CRAM -> VCF and more, orchestrates [GATK-SV WDL workflows](https://github.com/broadinstitute/gatk-sv) in Hail Batch.

## Seqr Loader

![uml](docs/seqr_loader.png)

Seqr Loader is a combination of the following [WARP](https://github.com/broadinstitute/warp)-inspired workflows:

* Single-sample whole genome and WES germline calling, FASTQ -> GVCF (BWA/DRAGMAP, GATK4 HaplotypeCaller)
* Single-sample calling with re-alignment, CRAM -> GVCF (we use [BAZAM](https://github.com/ssadedin/bazam) to extract FASTQ from CRAM/BAM)
* Whole genome and WES joint-calling, GVCFs -> pVCF
* WES AS-VQSR and WGS VQSR workflows

As well as:

* [the Broad's seqr loading pipelines](https://github.com/broadinstitute/seqr-loading-pipelines)
* [Somalier](https://github.com/brentp/somalier) pedigree checks
* [MultiQC](https://github.com/ewels/MultiQC) QC reporting.

The workflows use Metamist as a source of FASTQ, CRAMs, and sample/participant metadata, and TOML configs for extended configuration (dataset and samples, Hail Batch parameters, Elasticsearch credentials, QC thresholds).

## Installation

Requirements:

* [Analysis runner](https://github.com/populationgenomics/analysis-runner)

Clone the repository recursively and change into the repository folder:

```bash
git clone --recurse-submodules git@github.com:populationgenomics/production-pipelines.git
cd production-pipelines
```

## Usage

To run Seqr Loader, there is a script called `main.py`, which can take configs as input and should be submitted using the analysis-runner:

```bash
analysis-runner \
  --dataset seqr \
  --description "Seqr Loader" \
  --output-dir "seqr" \
  --access-level test \
  --config configs/genome.toml \
  --config configs/acute-care-test.toml
  main.py
```

## Cookbook

Seqr Loader can be used partially, controlled by `workflows/first_stages`, `workflows/last_stages`, or `workflows/only_stages` parameters.

[//]: # (You can use the to only do joint-calling &#40;GVCF -> pVCF&#41;, run AS-VQSR, run matrix table annotation, create Elasticsearch index, create MultiQC reports, etc., )

### Align Seqr genomes or exomes and produce CRAM QC

All configuration files that is passed with the `--config` option are combined by analysis-runner. In order to align data, you only need to trigger the `Align` stage from the workflow, which can be done with `workflows/first_stages`. So you can create a configuration file like the following:

```toml
[workflow]
last_stages = ['Align']
```

And assuming it's named `config.toml`, run:

```bash
analysis-runner --dataset seqr --description "Loader" --output-dir "loader" \
  --access-level test \
  --config configs/genome.toml \
  --config configs/seqr-main.toml \
  --config config.toml \
  main.py
```

For exomes, replace `configs/genome.toml` with `configs/exome.toml`, or set `sequencing_type = 'exome'` in the `workflow` section.

The section `workflows/input_datasets` in `configs/seqr-main.toml` specified the list of all projects to be processed, excluding samples specified in the `workflows/skip_samples` section.

QC reports for each dataset would be exposed on a web server, e.g. for `validation` genomes, the URL would be https://main-web.populationgenomics.org.au/validation/qc/cram/multiqc.html, and for exomes, it will  be https://main-web.populationgenomics.org.au/validation/exome/qc/cram/multiqc.html.

If samples had pedigree data, a Somalier report will be run to infer and validate participant relationships, with report produced as https://main-web.populationgenomics.org.au/validation/qc/cram/somalier.html

### Call GVCFs for each sample and validate variant calls

Do the same as above, but with the following section in `config.toml`:

```toml
[workflow]
first_stages = ['Genotype']
last_stages = ['Genotype']
```

The genome GVCF QC report will be exposed as https://main-web.populationgenomics.org.au/validation/qc/gvcf/multiqc.html, for `validation` samples it would include a [hap.py](https://github.com/Illumina/hap.py) section with validation stats.

### Upload Elasticsearch indices

If you want the workflow to create Elasticsearch indices in the end, run the entire workflow, but specify the `workflow/create_es_index_for_datasets` section with the list of datasets for which you want the indices to be created:

```toml
[workflow]
create_es_index_for_datasets = ['validation']
```

The resulting index will be named using the current datestamp, or using `worfklow/output_version` option if it's specified. The Elasticsearch server is configured using the `elasticsearch` section in `configs/seqr.toml`.

### Run the large cohort workflow

The `Genotype` stage can alternatively be followed by [the large cohort workflow](https://github.com/populationgenomics/large-cohort-pipeline) for:
	* VDS combiner to joint-call GVCFs directly into a VDS format,
	* Run Hail sample and variant QC, infer ancestry, relatedness, and sex,
	* Build CPA plots.

Refer to the `large-cohort-pipeline` repository for details.

## Development

### Dev installation

Requirements:

* python 3.10+

Clone the repository recursively and change into the repository folder:

```bash
git clone --recurse-submodules git@github.com:populationgenomics/production-pipelines.git
cd production-pipelines
```

Assuming `python3` points to `python3.10`, set up the environment:

```bash
python3 -m pip install virtualenv
virtualenv venv
source venv/bin/activate
pip install -r requirements-dev.txt
pip install -r requirements.txt
```

### Stages and jobs

The `stages` folder provides implementations of the Seqr-Loader stages, that can be used altogether or in isolation. Stages use functions in the `jobs` module, that create jobs for different applications, e.g. alignment, deduplication, fastqc, haplotype caller, creating calling intervals, splitting VCF, running VQSR, etc. Those jobs can be called from the `queue_jobs()` method of a stage, or in isolation, as they don't need to know about the stage they are called from, but only need a `Batch` object. Feel free to explore both `stages` and `jobs` folders for examples and inspiration.
