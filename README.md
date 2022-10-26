# Hail Batch Workflows

This repository provides implementations of genomics workflows using Hail Batch and the [cpg_utils/workflows library](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/workflows/README.md), specifically:

* WES/WGS Seqr Loader: FASTQ -> CRAM -> GVCF -> pVCF -> Matrix Table -> Elasticsearch index, with an ability to use parts of this pipeline as e.g. a single-sample germline variant calling workflow (FASTQ -> GVCF), joint-calling pipeline (GVCF -> pVCF), AS-VQSR pipeline, etc.
* WES/WGS Large Cohort Workflow: FASTQ -> CRAM -> GVCF -> VDS with Hail Table annotations.
* (in progress) GATK-SV: CRAM -> VCF and more, orchestrates [GATK-SV WDL workflows](https://github.com/broadinstitute/gatk-sv) in Hail Batch.

## Installation

Requirements:

* [Analysis runner](https://github.com/populationgenomics/analysis-runner)

Clone the repository recursively and change into the repository folder:

```bash
git clone --recurse-submodules git@github.com:populationgenomics/production-pipelines.git
cd production-pipelines
```

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

### Example usage

Run `main.py seqr_loader` through the analysis runner:

```bash
analysis-runner \
  --dataset seqr --description "Seqr Loader validation" --output-dir "seqr" \
  --access-level full \
  --config configs/genome.toml \
  --config configs/validation.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest \
  main.py seqr_loader
```

Note the configs that are being passed to analysis-runner: `configs/genome.toml` and `configs/validation.toml`. They are merged together by analysis-runner, with values in the latter overwriting values in the former. For more about configs, see [team-docs](https://github.com/populationgenomics/team-docs/blob/main/cpg_utils_config.md). Note that the `seqr_loader` command loads the [configs/defaults/seqr_loader.toml](configs/defaults/seqr_loader.toml) by default.

### Seqr production load invocation

`configs/seqr-main.toml` provides a configuration of a CPG production Seqr load: specifically, the list of datasets to process and joint-call together, and a list of blacklisted samples in those datasets. Another config, `configs/genome.toml` or `configs/exome.toml`, can be used to subset samples to WGS or WES specifically. One of these two must be provided, as the Seqr loader can work on only one type of data at a time. 

For example, to load the genome data:

```sh
analysis-runner \
  --dataset prophecy --description "Seqr Load" --output-dir "seqr" \
  --access-level full \
  --config configs/seqr-main.toml \
  --config configs/genome.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest \
  main.py seqr_load
```

#### Align Seqr genomes or exomes and produce CRAM QC

Seqr Loader can be used partially, controlled by `workflows/first_stages`, `workflows/last_stages`, or `workflows/only_stages` parameters.

In order to align data, you only need to trigger the `Align` stage from the workflow (or `CramMultiQC`, if you also want a dataset-level MultiQC report), which can be done with `workflows/last_stages`. So you can create a configuration file like the following:

```toml
[workflow]
last_stages = ['CramMultiQC']
```

And assuming it's named `~/myconfig.toml`, run:

```bash
analysis-runner \
  --dataset seqr --description "CRAM MultiQC" --output-dir "seqr" \
  --access-level full \
  --config configs/genome.toml \
  --config configs/validation.toml \
  --config ~/myconfig.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest \
  main.py seqr_loader
```

For exomes, replace `configs/genome.toml` with `configs/exome.toml`, or set `sequencing_type = 'exome'` in the `workflow` section.

The section `workflows/input_datasets` in `configs/seqr-main.toml` specified the list of all projects to be processed, excluding samples specified in the `workflows/skip_samples` section.

QC reports for each dataset would be exposed on a web server, e.g. for `validation` genomes, the URL would be https://main-web.populationgenomics.org.au/validation/qc/cram/multiqc.html, and for exomes, it will  be https://main-web.populationgenomics.org.au/validation/exome/qc/cram/multiqc.html.

If samples had pedigree data, a Somalier report will be run to infer and validate participant relationships, with report produced as https://main-web.populationgenomics.org.au/validation/qc/cram/somalier.html

### Call GVCFs for each sample and validate variant calls

Do the same as above, but with the following section in `~/myconfig.toml`:

```toml
[workflow]
last_stages = ['GvcfMultiQC']
```

The genome GVCF QC report will be exposed as https://main-web.populationgenomics.org.au/validation/qc/gvcf/multiqc.html, for `validation` samples it would include a [hap.py](https://github.com/Illumina/hap.py) section with validation stats.

### Upload Elasticsearch indices

If you want the workflow to create Elasticsearch indices in the end, run the entire workflow, but specify the `workflow/create_es_index_for_datasets` section with the list of datasets for which you want the indices to be created:

```toml
[workflow]
create_es_index_for_datasets = ['validation']
```

The resulting index will be named using the current datestamp, or using `worfklow/output_version` option if it's specified. The Elasticsearch server is configured using the `elasticsearch` section in `configs/seqr.toml`. The reason for not automatically creating indices for every project is that the Elasticsearch instance can easily run out of disk space, so additional safeguard is needed. 

## Large Cohort Workflow

![uml](docs/large_cohort.png)

A [Hail Query](https://hail.is/) workflow for large germline genomic variant calling cohorts.

1. Combine GVCFs (generated by GATK4) into a [VDS](https://hail.is/docs/0.2/vds/hail.vds.VariantDataset.html#hail.vds.VariantDataset) format.
2. Perform sample-level QC, including, pedigree, sex and ancestry inference.
3. Perform variant-level QC, including [allele-specific VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants).

### Usage

```sh
analysis-runner \
  --dataset prophecy --description "Larcoh thousand-genomes" --output-dir "larcoh" \
  --access-level test \
  --config configs/test.toml \
  --config configs/thousand-genomes.toml \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest \
  main.py large_cohort
```

The workflow will find GVCFs for input samples using Metamist, along with available sample metadata (e.g. known population labels, sex, QC), and would write the results into the `gs://cpg-prophecy-test` bucket.

Note that the `large_cohort` command loads the [configs/defaults/large_cohort.toml](configs/defaults/large_cohort.toml) by default, with other configs specified with `--config` put on top. For more about configs, see [team-docs](https://github.com/populationgenomics/team-docs/blob/main/cpg_utils_config.md). 

### Outputs

#### VDS

GVCFs are combined into a [VDS folder format](https://hail.is/docs/0.2/vds/hail.vds.VariantDataset.html#variantdataset) which is writen to `gs://cpg-prophecy-test/vds/v0-1.vds`.

#### Sample QC

* Sample-level metadata and QC is written to `gs://cpg-prophecy-test/large-cohort/v0-1/sample_qc.ht`, with the following row fields.

Metamist metadata:

```
's': str
'external_id': str
'dataset': str
'gvcf': str
'sex': str
'continental_pop': str
'subcontinental_pop': str
```

[hl.sample_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc) annotations:

```
'sample_qc': struct {
    n_het: int64,
    n_hom_var: int64,
    n_non_ref: int64,
    n_singleton: int64,
    n_singleton_ti: int64,
    n_singleton_tv: int64,
    n_snp: int64,
    n_insertion: int64,
    n_deletion: int64,
    n_transition: int64,
    n_transversion: int64,
    n_star: int64,
    r_ti_tv: float64,
    r_ti_tv_singleton: float64,
    r_het_hom_var: float64,
    r_insertion_deletion: float64,
    bases_over_gq_threshold: tuple (
        int64,
        int64,
        int64
    ),
    bases_over_dp_threshold: tuple (
        int64,
        int64,
        int64,
        int64,
        int64
    )
}
```

[Sex imputation](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.impute_sex):

```
'is_female': bool
'chr20_mean_dp': float64
'chrX_mean_dp': float64
'chrX_ploidy': float64
'chrY_mean_dp': float64
'chrY_ploidy': float64
'X_karyotype': str
'Y_karyotype': str
'sex_karyotype': str
'impute_sex_stats': struct {
    f_stat: float64,
    n_called: int64,
    expected_homs: float64,
    observed_homs: int64
}
```

Soft filters are populated based on the above results:

```
'filters': set<str>
```

Sample-QC based filters are calculated according to the thresholds specified in the config TOML, with available values `low_coverage` and `bad_sample_qc_metrics`:

```toml
[large_cohort.sample_qc_cutoffs]
min_coverage = 18
max_n_snps = 8000000
min_n_snps = 2400000
max_n_singletons = 800000
max_r_duplication = 0.3
max_r_het_hom = 3.3
```

Sex imputation based filters are `sex_aneuploidy` and `ambiguous_sex`.

#### Relatedness

[PC-Relate method](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) is used to identify pairs of the 1st and the 2nd degree relatives (kin coefficient threshold - below which samples are considered unrelated - is specified as `large-cohort.max_kin` in TOML). Pairwise sample relatedness matrix is written as a Hail table index by a tuple of sample IDs: `gs://cpg-prophecy-test/large-cohort/v0-1/relatedness.ht`

```
Row fields:
    'i': str
    'j': str
    'kin': float64
    'ibd0': float64
    'ibd1': float64
    'ibd2': float64
----------------------------------------
Key: ['i', 'j']
```

`gs://cpg-prophecy-test/large-cohort/v0-1/relateds_to_drop.ht` is a sample-level table which contains related samples to drop, with top ranking sample selected from each family. Sets of unrelated individuals are determined using Hail's [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html?highlight=maximal_independent_set#hail.methods.maximal_independent_set). 

```
Row fields:
    's': str
    'rank': int64
----------------------------------------
Key: ['s']
```

### Ancestry

PCA results are written into `gs://cpg-prophecy-test/large-cohort/v0-1/ancestry`:
  * `gs://cpg-prophecy-test/large-cohort/v0-1/ancestry/eigenvalues.ht`
  * `gs://cpg-prophecy-test/large-cohort/v0-1/ancestry/loadings.ht`
  * `gs://cpg-prophecy-test/large-cohort/v0-1/ancestry/scores.ht`
 
When there are samples with known `continental_pop` available, using the PCA results a random forest method is used to infer population labels. The method is trained using 16 principal components as features on samples with known ancestry. Ancestry was assigned to all samples for which the probability of that ancestry was high enough (the threshold is configured as `large_cohort.min_pop_prob` in TOML). Results are written as sample-level table `gs://cpg-prophecy-test/large-cohort/v0-1/ancestry/inferred_pop.ht`.

```
Row fields:
    's': str
    'scores': array<float64>
    'pop': str
    'is_training': bool
    'pca_scores': array<float64>
----------------------------------------
Key: ['s']
``` 

Plots for PCA and loadings are written to `gs://cpg-prophecy-test-web/large-cohort/v0-1/ancestry/*`, a bucket that is exposed as https://test-web.populationgenomics.org.au/prophecy/large-cohort/ancestry/*
 
#### Dense subset

For PCA and PC-relate, a dense subset of the original dataset is used. The markers for the subset are read from the `references.gnomad.predetermined_qc_variants` Hail table specified in the TOML config. The table is suitable for both exomes and genomes, so that a mixture of different sequencing types can be processed together. The resulting subset is written to `gs://cpg-prophecy-test/large-cohort/v0-1/dense-subset.mt`.

### Allele-specific variant quality score recalibration (AS-VQSR)

Variants from good quality samples are filtered using the [AS-VQSR method](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-):

1. Variants are exported into a sites-only VCF,

1. Batch jobs are submitted to create SNP and indel recalibration models using the allele-specific version of GATK Variant Quality Score Recalibration [VQSR](https://gatkforums.broadinstitute.org/gatk/discussion/9622/allele-specific-annotation-and-filtering), with the standard GATK training resources (HapMap, Omni, 1000 Genomes, Mills indels), and the following features:
   
   * SNVs:   `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`, `AS_MQ`,
   * Indels: `AS_FS`, `AS_SOR`, `AS_ReadPosRankSum`, `AS_MQRankSum`, `AS_QD`.
   
1. The models are applied to the VCFs and combine them back into one VCF.
   
1. VCF is converted back into a sites-only locus-level Hail table `gs://cpg-prophecy-test/large-cohort/v0-1/vqsr.ht`, with split multiallelics.

```
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'filters': set<str>
    'info': struct {
        NEGATIVE_TRAIN_SITE: bool,
        POSITIVE_TRAIN_SITE: bool,
        culprit: str
    }
    'a_index': int32
    'was_split': bool
```

Note that the `info.AS-*` annotations used for AS-VQSR are dropped, and only the resulting filter label is appended into the `filters` field, e.g. `VQSRTrancheINDEL99.50to99.90`, `VQSRTrancheSNP99.00to99.90+`, etc. The AS_VQSLOD thresholds for assigning filters are configurable in TOML as `vqsr.snp_filter_level` and `vqsr.indel_filter_level`.

This pipeline is largely compiled from the following two WDL workflows:
   
1. `hail-ukbb-200k-callset/GenotypeAndFilter.AS.wdl`

2. The [Broad VQSR workflow](https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl) documented [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering), translated from WDL with a help of [Janis](https://github.com/PMCC-BioinformaticsCore/janis).

#### Frequencies

Frequencies are calculated using the Hail's [hl.variant_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.variant_qc) method from good quality samples, and written to `gs://cpg-prophecy-test/large-cohort/v0-1/frequencies.ht` locus-level table with split multiallelics:

```
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'a_index': int32
    'was_split': bool
    'InbreedingCoeff': float64
    'dp_stats': struct {
        mean: float64,
        stdev: float64,
        min: float64,
        max: float64
    }
    'gq_stats': struct {
        mean: float64,
        stdev: float64,
        min: float64,
        max: float64
    }
    'AC': array<int32>
    'AF': array<float64>
    'AN': int32
    'homozygote_count': array<int32>
    'call_rate': float64
    'n_called': int64
    'n_not_called': int64
    'n_filtered': int64
    'n_het': int64
    'n_non_ref': int64
    'het_freq_hwe': float64
    'p_value_hwe': float64
    'p_value_excess_het': float64
----------------------------------------
Key: ['locus', 'alleles']
```

### Applied outputs

Tables generated by the workflow can be applied to a VDS or a split+dense matrix table in the following way:

```python
import hail as hl

vds = hl.vds.read_vds('gs://cpg-prophecy-test/vds/v0-1.vds')
sample_qc_ht = hl.read_table('gs://cpg-prophecy-test/large-cohort/v0-1/sample_qc.ht')
relateds_to_drop_ht = hl.read_table('gs://cpg-prophecy-test/large-cohort/v0-1/relateds_to_drop.ht')
pop_ht = hl.read_table('gs://cpg-prophecy-test/large-cohort/v0-1/ancestry/inferred_pop.ht')
vqsr_ht = hl.read_table('gs://cpg-prophecy-test/large-cohort/v0-1/vqsr.ht')
freq_ht = hl.read_table('gs://cpg-prophecy-test/large-cohort/v0-1/frequencies.ht')

# Row-level tables require a split+dense matrix table:
vds = hl.vds.split_multi(vds, filter_changed_loci=True)
mt = hl.vds.to_dense_mt(vds)

# Hard-filtering samples and variants:
mt = mt.filter_cols(hl.len(sample_qc_ht[mt.col_key].filters) > 0, keep=False)
mt = mt.filter_cols(hl.is_defined(relateds_to_drop_ht[mt.col_key]), keep=False)
mt = mt.filter_rows(hl.len(vqsr_ht[mt.row_key].filters) > 0, keep=False)

# Annotating samples and variants:
mt = mt.annotate_cols(**sample_qc_ht[mt.col_key])
mt = mt.annotate_cols(**pop_ht[mt.col_key])
mt = mt.annotate_rows(**vqsr_ht[mt.row_key])
mt = mt.annotate_rows(**freq_ht[mt.row_key])
```

### References

The workflow is largely inspired by [the Hail pipeline used for the QC of gnomAD releases](https://github.com/broadinstitute/gnomad_qc). Good summaries of gnomAD QC can be found in gnomAD update blog posts:

* [https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad](https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad)
* [https://macarthurlab.org/2018/10/17/gnomad-v2-1](https://macarthurlab.org/2018/10/17/gnomad-v2-1)
* [https://macarthurlab.org/2019/10/16/gnomad-v3-0](https://macarthurlab.org/2019/10/16/gnomad-v3-0)
* [https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control](https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#sample-and-variant-quality-control)
* [https://blog.hail.is/whole-exome-and-whole-genome-sequencing-recommendations/](https://blog.hail.is/whole-exome-and-whole-genome-sequencing-recommendations/)

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

### `cpg-utils/workflows`

The codebase uses the `cpg-utils/workflows` library to drive workflows. If you are making changes to that library, you can create a pull request for [cpg-utils](https://github.com/populationgenomics/cpg-utils), which would trigger a Docker build, tagged with a Git commit SHA. Then you can pass that image to the analysis runner with `--image`, in order to use your `cpg-utils` change without pushing it to main.

```bash
analysis-runner --dataset thousand-genomes \
--description "thousand-genomes larcoh" --output-dir "thousand-genomes" --access-level test \
--config configs/genome.toml --config configs/thousand-genomes.toml \
--image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:b4cadbc63f406345c67437fd10efb3a6108f1726 \
main.py large_cohort
````
