# [Large Cohort Pipeline](https://github.com/populationgenomics/production-pipelines)
## Quick Overview
This pipeline takes in FASTQ or CRAM files as input and runs a series of python scripts to process all files into a final annotated, joint-called [sparse matrix table](https://hail.is/docs/0.2/vds/index.html ""). Specifically, the pipeline takes on the following stages:
1. CramMultiQC Stage:  Aligns the data and generates a MultiQC report
2. GvcfMultiQC Stage: Generates GVCFs and generates a MultiQC report
3. FrequenciesStage: Runs the [VDS combiner](https://hail.is/docs/0.2/vds/hail.vds.combiner.VariantDatasetCombiner.html#hail.vds.combiner.VariantDatasetCombiner "") with GVCFs as an input to produce a VariantDataset (VDS), then performs sample-level QC (pedigree, sex and ancestry inference)
4. LoadVqsr Stage: Converts the VDS into a sites-only VCF-ready table, then performs allele-specific VQSR variant QC.

Graphically, the pipeline has the following workflow:

![uml](large_cohort.png)

## Running the pipeline
To run, the pipeline takes in at least one config file, however multiple can be given. As an example, an exome/genome config file can be provided specifying whether samples are whole genome or whole exome sequences (e.g., `configs/genome.toml`), and a config file specifying parameters specific to the dataset can be provided with any name the user chooses. Any additional config files can be provided by adding  `--config` to the analysis runner command. All config files are then merged with a default config file,  `configs/defaults/large_cohort.toml`, to form one master file which is fed into the analysis. Any duplicated entries in separate config files will result in entries from the last config file overwriting previous ones. 

An example run would look like the following:

```sh
analysis-runner --dataset bioheart --description "bioheart larcoh" \
--output-dir "bioheart" --access-level test \
--config configs/genome.toml --config configs/bioheart-test.toml \
--image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest main.py large_cohort
```

## Detailed stage descriptions
The detailed description for each stage is as follows:
### CramMultiQC Stage

1.  Alignment: the first step is to align the given input files. By default, DRAGMAP is used, with the version matching the tag of the corresponding image used for the alignment job (specified in. [images] within the config file). This then outputs a CRAM file with a corresponding index file (.crai file). 
<ins>Configurable inputs:</ins>
	* [workflow].check_intermediates
	* [workflow].skip_samples_with_missing_input
	* [workflow].cram_version_reference
	* [workflow].realign_from_cram_version
	* [workflow].check_inputs
	* [workflow].sequencing_type
	* [workflow].resources.Align.storgae_gb

2. CramQC: This step calls a set of tools that process CRAM files for QC purposes. It relies on the jobs `somalier.py`, `picard.py`, `samtools.py`, and `verifybamid.py` to run qc metrics on all samples and is based off of the [Broad Warp pipeline](https://github.com/broadinstitute/warp/tree/master ""). Specifically, the following QC functions are run:
	1. `somalier.extract`: this generates a fingerprint using [Somalier](https://github.com/brentp/somalier ""), with the output being a `*.somalier` file.
	2.  `picard_wgs_metrics`: This runs picard CollectWgsMetrics metrics, which collects metrics about the fractions of reads that pass base and mapping-quality filters, as well as coverage (read-depth) levels for WGS analyses (see the [CollectWgsMetrics help page](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard- "")). 
	3.  `picard_hs_metrics`: This take a SAM/BAM file as input and collects Hybrid-selection (HS) metrics, as described in the [Picard metric definitions help page.](http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics "") 
	4.  `picard_collect_metrics`: This runs Picard’s CollectMultipleMetrics, which is a ‘meta-metrics' tool that runs one or more metrics collection modules at the same time to cut down on the time spent reading in data from input files. For all available modules, see the [Picard CollectMultipleMetrics help page](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard- ""). The tool produces outputs of '.pdf' and '.txt' files for each module, except for the `CollectAlignmentSummaryMetrics module`, which outputs only a '.txt' file.
	5.  `samtools_stats`: this runs  `samtools stats` for alignment QC, which produces a comprehensive list of statistics from the alignment file. For all output statistics, see the [samtools-stats page](http://www.htslib.org/doc/samtools-stats.html "")
	6.  `verifybamid`: this runs `VerifyBamID` contamination checks, as explained in the [VerifyBamID GitHub page](https://github.com/Griffan/VerifyBamID ""). The expected output is a *.selfSM file, along with a TSV file with 2 rows and 19 columns. The first row contains keys (e.g., `SEQ_SM`, `RG`, `FREEMIX`) while the second row contains associated `VerifyBamID` values.
<br><ins>Configurable inputs:</ins>
	* [workflow].dry_run
	* [workflow].sequencing_type
	* [workflow].check_intermediates

3. SomalierPedigree: This step runs the `somalier.pedigree` function, which adds somalier and peddy jobs that infer relatedness and sex. These are then compared to the provided PED file with an attempt to recover them. If they are unable to be recoverd, further workflow jobs are cancelled. This step relies on the job `somalier.py`, with an output report being generated by [MultiQC](https://multiqc.info/ ""). 
<ins>Configurable inputs:</ins>
	* NA

4. CramMultiqc: This step runs [MultiQC](https://multiqc.info/ "") to aggregate CRAM QC stats (from above), using the `multiqc.py` job to generate the output.
<br><ins>Configurable inputs:</ins>
	* [workflow].skip_qc
	* [qc_thresholds]

### GvcfMultiQC Stage
1. Genotype: This step uses the `genotype` job to run the HaplotypeCaller, which genotypes individual samples (i.e. CRAM -> GVCF) and produces a GVCF + corresponding TBI index. 
<br><ins>Configurable inputs:</ins>
	* [workflow].check_intermediates
	* [workflow].intervals_path
	* [workflow].[sequencing_type]

2. GvcfQC: This step calls tools that process GVCFs for QC purposes. It runs Picard’s `CollectVariantCallingMetrics`, which collects metrics relating to snps and indels within a variant-calling file (see the [CollectVariantCallingMetrics help page](https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectVariantCallingMetrics.VariantCallingDetailMetrics "")). This step is reliant on the previous `genotype` step, but run at a parallel stage with `GvcfHappy`. 
<br><ins>Configurable inputs:</ins>
	* [workflow].skip_qc
	* [workflow].check_intermediates
	* [workflow].[sequencing_type]

3. GvcfHappy: This step runs [Hap.py](https://github.com/Illumina/hap.py/blob/master/doc/happy.md "") validation/concordance stats to validate a GVCF for samples where a truth callset is available. The expected output is a MultiQC `*.summary.csv` file. 
<br><ins>Configurable inputs:</ins>
	* [workflow].[sequencing_type]
	* [validation].sample_map

4. GvcfMultiQC: This step runs MultiQC to summarise all GVCF QC output, with an HTML and corresponding JSON file produced.
<ins>Configurable inputs:</ins>
	* [workflow].skip_qc
	* [workflow].dry_run
	* [qc_thresholds]
### FrequenciesStage
1. Genotype: This step uses the `genotype` job to run the HaplotypeCaller, which genotypes individual samples (i.e. CRAM -> GVCF) and produces a GVCF + corresponding TBI index. 
<br><ins>Configurable inputs:</ins>
	* [workflow].check_intermediates
	* [workflow].intervals_path
	* [workflow].[sequencing_type]

2. Combiner: This step runs the Hail [VDS combiner](https://hail.is/docs/0.2/vds/hail.vds.combiner.VariantDatasetCombiner.html#hail.vds.combiner.VariantDatasetCombiner ""), which produces a VariantDataset (VDS) by combining any number of GVCF and/or VDS files. For more information, see Hail’s documenatation on [Variant Datasets](https://hail.is/docs/0.2/vds/index.html ""). The output of this step is written to `gs://cpg-{dataset}-{access-level}/vds/v0-1.vds`.
<br><ins>Configurable inputs:</ins>
	* [workflow].vds_version
	* [hail].dataproc.combiner_autoscaling_policy
	* get_config().get('combiner', {})
	* [workflow].[sequencing_type]
	* [workflow].check_inputs
	* [workflow].skip_samples_with_missing_input

3. SampleQC: This step takes an input VDS table and performs the following steps: 1. removes centromeres and telomeres 2. filters to autosomes and 3. computes sample quality metrics about the VDS (supplied in the `large_cohort.sample_qc_cutoffs` within the `large_cohort.toml` file). Sex is then imputed on the VDS and the following soft filters added: samples with ambiguous sex assignments, low-coverage samples, and extreme raw bi-allelic sample QC outliers. Sample-QC based filters are calculated according to the thresholds specified in the `large_cohort.sample_qc_cutoffs` section of the `configs/defaults/large_cohort.toml`. The output of this step is written to `gs://cpg-{dataset}-{access-level}/large-cohort/v0-1/sample_qc.ht`. 
<ins>Configurable inputs:</ins>
	* [workflow].[sequencing_type]
	* [workflow].[sequencing_type]
	* [large_cohort].[sample_qc_cutoffs]

4. DenseSubset: This step filters a sparse VDS to a set of predetermined QC sites and returns a dense (no filtered entries) MatrixTable with split multiallelics. The markers for the subset are read from the `references.gnomad.predetermined_qc_variants` Hail table within the `configs/defaults/large_cohort.toml` config. The table is suitable for both exomes and genomes, so that a mixture of different sequencing types can be processed together. The resulting subset is written to `gs://cpg-{dataset}-{access-level}/large-cohort/v0-1/dense-subset.mt`.
<br><ins>Configurable inputs:</ins>
	* NA

5. Relatedness: This step takes a hail matrix table as input and runs [pcrelate](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4716688/ ""). Under the PC-Relate model, kinship ranges from 0 to 0.5, where 0.5 = identical (i.e., identical twins, or the same individual), parent-child and sibling pairs have a kinship of 0.25, and uncle_aunt_grand-parent/-child pairs have a kinship of 0.125. Related individuals in the config file are defined as pairs of 1st and 2nd degree relatives, which can be adjusted in the `large-cohort.max_kin` argument in the `large_corhort.toml` file. This will write a table with the following structure:

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
A pairwise sample relatedness matrix is written as a Hail table to: `gs://{project-name}/relatedness.ht`.
Once pcrelate is performed, a sample-level table which contains related samples to drop is created using Hail’s [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html?highlight=maximal_independent_set#hail.methods.maximal_independent_set)and stored in `gs://{project-name}/relateds_to_drop.ht`.
<br><ins>Configurable inputs:</ins><br>
* [large_cohort].[max_kin]

6. Ancestry: This runs a PCA, excluding given related samples, and projects samples in PC space to return their scores (taken from the `run_pca_with_relateds` [function from gnomAD](https://broadinstitute.github.io/gnomad_methods/api_reference/sample_qc/ancestry.html#gnomad.sample_qc.ancestry.run_pca_with_relateds "")). This function returns four tables: `scores.ht`, `eigenvalues.ht`, `loadings.ht`, and an original sample ht (annotated with the `training_pop` , `pca_scores`, `pop`, and `prob` for each population label). The output is stored in `gs://{project-name}/ancestry` bucket. 
When there are samples with known `continental_pop` available, a random forest method is used to infer population labels. The method is trained using 16 principal components as features on samples with known ancestry. Ancestry is assigned to all samples for which the probability of that ancestry is high enough (the threshold is configured as `large_cohort.min_pop_prob` in TOML). Results are written as sample-level table `gs://cpg-{dataset}-{access-level}/large-cohort/v0-1/ancestry`.
<br><ins>Configurable inputs:</ins>
	* [large_cohort].[min_pop_prob]
	* [large_cohort].[n_pcs]

7. AncestryPlots: This plots a scatter plot of the ancestry scores and loadings. The output is an HTML file, written for each PC and scope ("dataset", "population”).
<br><ins>Configurable inputs:</ins>
	* NA

8. Frequencies: Generates frequency annotations (AF, AC, AN, InbreedingCoeff) for an input VDS using Hail’s [hl.variant_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.variant_qc) function. The output is written with split multiallelics to `gs://cpg-{dataset}-{access-level}/large-cohort/v0-1/frequencies.ht`.
<br><ins>Configurable inputs:</ins>
	* NA

### LoadVqsr Stage
1. MakeSiteOnlyVcf: This step converts a VDS into a sites-only VCF-ready table.
<br><ins>Configurable inputs:</ins>
	* NA

2. Vqsr: This step adds jobs that perform the allele-specific VQSR variant QC
<br><ins>Configurable inputs:</ins>
	* [workflow].use_as_vqsr
	* [workflow].get(‘intervals_path’)
	* [vqsr].[snp_filter_level]
	* [vqsr].[indel_filter_level]

3. LoadVqsr: This ste converts VQSR VCFs to HTs.
<br><ins>Configurable inputs:</ins>
	* NA

## Notes on what we need
* all parameters should have options: e.g., Options available are `runner`, `text`, and `json`
* All thresholds made should have justification as to why that was chosen
	* Curious as to why `min_pop_prob` was chosen to be 0.5. In gnomad v3.0, this is set to 90% (see https://macarthurlab.org/2019/10/16/gnomad-v3-0/)
	* Why is 0.2 used as the unrelated threshold? Avuncular pairs and grand-parent/child pairs both have kinship 0.125 in expectation and both have identity-by-descent-zero 0.5 in expectation
* Any unused values from config files should be removed
	* scatter_count
	* max_r_duplication
	* gencode_gtf = 'gencode/gencode.v39.annotation.gtf.bgz'
	* references.gnomad section
* How is get_config()['workflow'] differentiated between the bioheart-test.toml and the large_cohort.toml?
	* I don’t actually know which toml file arguments are being pulled from, e.g., what differentiates `workflow.sequencing_type` from the `bioheart-test.toml` and `workflow.sequencing_type` from the `genome.toml`?
* How do you set up the number of workers, cores, etc in a run? Where does `scatter_count` feed into the pipeline?
* Why does `[workflow].sequencing_type` need to be specified in a separate exome/genome config file? Can this not be added in one of the other toml files?
* Could the VQSR section in the README and `large_cohort.toml` be described better? I have no idea what these are doing, and why these inputs are chosen. Also, is this only to be specified in the `large_cohort.toml`?
* Where is `get_config().get()` being pulled from?

## To do:
- [ ] ensure all correct default values are added in the configurable inputs parameters field
- [ ] make all descriptions in configurable inputs easy to understand (some of these are still hard to follow)
- [x] redistribute settings in `unsure` configurable inputs section to corresponding toml files

## Configurable inputs

All inputs below are either required or optional configurable inputs to the large cohort pipeline. Values can be entered into one or multiple config files, which are then merged with the `large_cohort.toml` found in the `configs/defaults/` folder. If an entry in one of the config files provided by the user is duplicated, the entry from the last config file will get used. E.g.,  if [large_cohort.sample_qc_cutoffs].min_coverage within the `large_cohort.toml` is set to 18 but set to 20 within a line of a later config supplied by the user, 20 will be used. 

| Setting                                           | Parameters                                                   | Description                                                  |
|---------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|
| [workflow].check_intermediates                    | boolean value; default = true                                | Within jobs, check all in-job intermediate files for possible reuse. If set to False, this will overwrite all intermediates. |
| [workflow].skip_qc                                | boolean value; default = false                               | Whether or not to skip the QC stage of the pipeline (==specifically which QC stage??)== |
| [workflow].use_as_vqsr                            | boolean value; default = true                                | Use allele-specific annotations for VQSR<br>use_as_vqsr = true — ==could this be described in a way that’s easy for the user to understand? I don’t fully understand what’s happening here== |
| [workflow].vds_version                            | Whole number or a floating point number (e.g., ‘0.1’)        | A numerical version to assign as the version number for the analysis run.  |
| [workflow].cram_version_reference                 | ==what are the available options?==                          | ==Need info==                                                |
| [workflow].check_inputs                           | boolean value; default = true                                | Check input file existence (e.g. FASTQ files). If files are missing, the --skip-samples-with-missing-input option controls whether such samples should be ignored, or raise an error |
| [workflow].intervals_path                         | defauls to whole genome intervals; ==what are the other options?== | Calling intervals (defauls to whole genome intervals)        |
| [workflow].dry_run                                | ==What is the default, if there is one?==                    | Only print the final merged config and a list of stages to be submitted.<br>This will skip any communication with Metamist, Hail Batch, and Cloud Storage, so the code can be run without permissions. |
| [workflow].input_datasets                         | default = NA                                                 | The dataset used for the large cohort pipeline run. This feeds into `inputs.py` |
| [workflow].realign_from_cram_version              | v + version number (integer), e.g., v1. The default value is v0. | Realign CRAM when available, instead of using FASTQ. The parameter value should correspond to CRAM version (e.g. v0 in `gs://cpg-fewgenomes-main/cram/v0/CPG01234.cram`) |
| [workflow].skip_samples_with_missing_input        | boolean value; default = true                                | For the first (not-skipped) stage, if the input for a target does not exist, just skip this target instead of failing. E.g. if the first stage is Align, and `sample.alignment_input` for a sample does not exist, remove this sample instead of failing. In other words, ignore samples that are missing results from skipped stages. |
| [workflow].input_datasets                         |                                                              |                                                              |
| [workflow].sequencing_type                        | choose between ‘genome’ or ‘exome’ as input                  | Specifies whether samples are whole genome or whole exome sequences |
| [workflow].resources.Align.storage_gb             | ==what are the available options?==                          | ==Need info - what is this doing?==                          |
| [somalier].exclude_high_contamination             | ==what are the available options?==                          | ==Need info==                                                |
| [validation].sample_map                           | ==Is there a default?==                                      | Map internally used validation sample external_id to truth sample names. ==This needs to be explained better here== |
| [vqsr].snp_filter_level                           | ==what are the available options? why is 99.7 chosen as the default?== | ==Need info==                                                |
| [vqsr].[indel_filter_level]                       | ==what are the available options? why is 99.0 chosen as the default?== | ==Need info==                                                |
| get_config().get('combiner', {})                  | ==what are the available options?==                          | ==What is this doing, and where is it coming from in the config?== |
| [qc_thresholds]                                   | ==what are the available options?==                          | additional QC thresholds provided by the user ==Need info==  |
| [large_cohort].n_pcs                              | numeric values of 1 - n features; default = 16               | Number of principal components to use for PCA in ancestry analysis (job ancestry_pca.py) |
| [large_cohort].min_pop_prob                       | Float value from 0.0 - 1.0; default = 0.5                    | The minimum random forest probability to use for population assignment |
| [large_cohort].max_kin                            | A float value from 0 - 0.5, where 0 is x and 0.5 is identical (monozygotic twins). The default value to be considered unrelated is 0.2 | The maximum kin threshold, caculated by [PC relate](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4716688/ ""), to be considered unrelated. Under the PC-Relate model, kinship, ϕijϕij​, ranges from 0 to 0.5, and is precisely half of the fraction-of-genetic-material-shared. Information on relatedness and relatedness thresholds can be viewed in the hail docs [here](https://hail.is/docs/0.2/methods/relatedness.html ""). |
| [large_cohort].pop_meta_field                     | ==what are the available options?==                          | ==Needs a description here==                                 |
| [large_cohort.sample_qc_cutoffs].min_coverage     | 18                                                           | ==Why was this threshold used? What is this used for?==      |
| [large_cohort.sample_qc_cutoffs].max_n_snps       | 8000000                                                      | ‘’                                                           |
| [large_cohort.sample_qc_cutoffs].min_n_snps       | 2400000                                                      | ‘’                                                           |
| [large_cohort.sample_qc_cutoffs].max_n_singletons | 800000                                                       | ‘’                                                           |
| [large_cohort.sample_qc_cutoffs].max_r_het_hom    | 3.3                                                          | ‘’                                                           |
| [hail].pool_label                                 | default = ‘large-cohort’                                     | no idea what this does                                       |
| [hail].combiner_autoscaling_policy                | ==xx in the format …==; default = vcf-combiner-50            | Autoscaling policy must be created in the project that corresponds to the analysis dataset. — ==this needs more info, and to be easier to understand== |
| [hail].delete_scratch_on_exit                     | boolean value; default = false                               | ==this needs more info==                                     |
