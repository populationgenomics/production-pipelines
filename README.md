# Production pipelines

## cpg_pipes

The `cpg_pipes` package provides Python functions that help to design a pipeline powered by [Hail Batch](https://hail.is/docs/batch/service.html), in the context of the CPG infrastructe:
- Google Cloud,
- [CPG storage policies](https://github.com/populationgenomics/team-docs/tree/main/storage_policies),
- The [analysis-runner](https://github.com/populationgenomics/analysis-runner) as a proxy for submissions
- The [sample-metadata](https://github.com/populationgenomics/sample-metadata) database as a source of input data.

### Installation

Requires Python 3.10.

```bash
pip install cpg_pipes
```

### Building pipelines

A common use case of a bioinformatics pipeline at the CPG is a set of stages, with intermediate outputs stored on buckets according to the storage policies, and optionally tracked by the sample-metadata database.

The package implements a concept of _stage_, that can act on a _target_: e.g. a _sample_, a _project_, or a _cohort_. For example, a stage that runs whole genome alignment to produce a CRAM file, would act on a sample, and a stage that performs joint-calling would act on an entire cohort. 

Each stage declares paths to the outputs it would produce by overriding abstract `expected_result()` method; and defines how jobs are added into Hail Batch using `queue_jobs()` method. 

To declare a stage, inherit a class from SampleStage, ProjectStage, or CohortStage, implement the abstract methods, and wrap the class with a `@stage` decorator:

```python
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput
from cpg_pipes.pipeline.sample import Sample

@stage
class WriteSampleName(SampleStage):
    def expected_result(self, sample: Sample):
        return 'gs://cpg-thousand-genomes-test/cram/NA12878.cram'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        expected_path = self.expected_result(self.pipe)
        # This Job just writes the sample name to a file:
        j = self.pipe.b.new_job('Write sample name')
        j.command(f'echo {sample.id} > {j.output}')
        self.pipe.b.write_output(j.output, expected_path)
        # Construct StageOutput object:
        return self.make_outputs(sample, data=expected_path, jobs=[j])
```

The `queue_jobs` method is expected to retuns output of type `StageOutput`: you can call `self.make_outputs()` to construct that object.

Stages can depend on each other. Use the `requires_stages` parameter to `@stage` to set dependencies, and use the `inputs` parameter in `queue_jobs` to get the output of the previous stage:

```python
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput
from cpg_pipes.pipeline.sample import Sample

@stage(requires_stages=[WriteSampleName])
class ReadSampleName(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        sample_name_path = inputs.as_path(sample, stage=WriteSampleName)
        <...>
```

Stages can also communicate with the sample-metadata database to read or write outputs. E.g. a stage that calls a HaplotypeCaller could write results as Analysis entries of type "gvcf":

```python
from cpg_pipes.smdb.types import AnalysisType
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.stage import SampleStage

@stage(sm_analysis_type=AnalysisType.GVCF)
class HaplotypeCaller(SampleStage):
    <...>
```

Stage of differnet levels can depend on each other, and `cpg_pipes` will resolve that correctly. E.g. joint-calling taking GVCF outputs to produce a cohort-level VCF:

```python
from cpg_pipes.smdb.types import AnalysisType
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.stage import SampleStage, CohortStage, StageInput, StageOutput
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.cohort import Cohort

@stage(sm_analysis_type=AnalysisType.GVCF)
class HaplotypeCaller(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        job = add_haplotype_calling_job()
        return self.make_outputs(sample, data=self.expected_result(sample), jobs=[job])

@stage(sm_analysis_type=AnalysisType.JOINT_CALLING, requires_stages=HaplotypeCaller)
class JointCalling(CohortStage):
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        # Get outputs from previous stage. Because the HaplotypeCaller stage
        # acts on a sample, we use a method `as_path_by_target` that returns
        # a dictionary index by sample ID:
        gvcf_by_sid = inputs.as_path_by_target(stage=HaplotypeCaller)
        job = add_joint_calling_job(gvcf_by_sid)
        return self.make_outputs(cohort, data=self.expected_result(cohort), jobs=[job])
```

To submit the constructed pipeline to Hail Batch, initialise the `Pipeline` object, and call `submit_batch()`. You need to pass a name, a title, a version of a run, the namespace according to the storage policies (test/main/tmp), and the `analysis_project` name that would be used to communicate with the sample-metadata DB.

```python
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.namespace import Namespace

pipeline = Pipeline(
    name='my_pipeline',
    title='My Pipeline',
    output_version='v0-1',
    namespace=Namespace.TEST,
    analysis_project='fewgenomes',
)
pipeline.submit_batch()
```

### Sample-metadata database

Stage targets are organised in a form of projects and samples, folowing the model used in the sample-metadata database. Additional target _cohort_ is used to define all samples from all projects. The Pipeline class can populate correspondning targets if you provide input project names with `input_projects`:

```python
pipeline = Pipeline(
    <...>,
    input_projects=['hgdp', 'thousand-genomes'],
)
```

Under the hood, it will create a `cohort` field of type `Cohort`, which incapsulates a list of projects of type `Project`, each of which contains a list of samples of type `Sample`:

```
>>> pipeline.cohort.get_projects()
['perth-neuro', 'acute-care']

>>> for project in pipeline.cohort.get_projects():
>>>     samples = project.get_samples()
>>>     print(samples[0])
Sample(id='CPG68197', external_id='HGDP00001', project=hgdp, meta={'project': 'hgdp'}, alignment_input=AlignmentInput(fqs1=None, fqs2=None, bam_or_cram_path='gs://cpg-nagim-test/cram/HGDP00001.cram', index_path='gs://cpg-nagim-test/cram/HGDP00001.cram.crai'), pedigree=None),
...
```

If `Participant` entry is available, `sample.participant_id` will be populated. If corresponding `Sequence` is available, `sample.seq` will be populated, and `reads` metadata will be parsed and populated as `sample.alignment_input`. If corresponding `Analysis` entries exist, they wil be populated as `sample.sanalysis_by_type`. If `Family` data is available, it will be parsed and populated as `sample.pedigre`.

Communication with the sample-metadata DB is organised through the `cpg_pipes.smdb` module, a wrapper around the sample-metadata database API. E.g., to get a dictionary of Sample entries indexed by project name, run:

```python
from cpg_pipes.smdb import SMDB
from cpg_pipes.namespace import Namespace

smdb = SMDB(analysis_project='seqr')

samples_by_project = smdb.get_samples_by_project(
    project_names=['acute-care', 'perth-neuro'],
    namespace=Namespace.TEST,
)
```

To add an Analysis entry, run:

```python
from cpg_pipes.smdb import SMDB

smdb = SMDB(analysis_project='seqr')

analysis_id = smdb.create_analysis(
    type_='cram',
    output='gs://cpg-acute-care-test/cram/CPG12345.cram',
    status='completed',
    sample_ids=['CPG12345'],
    project_name='acute-care',
)
```

To add a Batch job that updates an Analysis entry, run:

```python
from cpg_pipes.smdb import SMDB

smdb = SMDB(analysis_project='seqr')

j = smdb.make_sm_completed_job(
    batch,
    analysis_id=analysis_id,
    analysis_type='cram',
    project_name='seqr',
    sample_name=['CPG12345'],
)
```

### Jobs

The `cpg_pipes.jobs` module defines functions that create Hail Batch Jobs for differrent bioinformatics purposes: alignment, fastqc, deduplication, variant calling, VQSR, etc. E.g. to implement the joint calling stage above, you can use:

```python
from cpg_pipes.smdb.types import AnalysisType
from cpg_pipes.pipeline.pipeline import stage
from cpg_pipes.pipeline.stage import SampleStage, CohortStage, StageInput, StageOutput
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.cohort import Cohort
from cpg_pipes.jobs import haplotype_caller, joint_genotyping

@stage(sm_analysis_type=AnalysisType.GVCF)
class HaplotypeCaller(SampleStage):
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput: 
        cram_path = inputs.as_path(target=sample, stage=Align)
        expected_path = self.expected_result(sample)
        job = haplotype_caller.produce_gvcf(
            b=self.pipe.b,
            output_path=expected_path,
            cram_path=cram_path,
            number_of_intervals=50,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[job])

@stage(sm_analysis_type=AnalysisType.JOINT_CALLING, requires_stages=HaplotypeCaller)
class JointCalling(CohortStage):
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        gvcf_by_sid = inputs.as_path_by_target(stage=HaplotypeCaller)

        expected_path = self.expected_result(cohort)
        job = joint_genotyping.make_joint_genotyping_jobs(
            b=self.pipe.b,
            out_vcf_path=expected_path,
            samples=cohort.get_all_samples(),
            gvcf_by_sid=gvcf_by_sid,
            depends_on=inputs.get_jobs(),
            smdb=self.pipe.get_db(),
        )
        return self.make_outputs(cohort, data=expected_path, jobs=[job])
```

Available jobs include alignment:

```python
from cpg_pipes.jobs import align
j align.align(
    b=b,
    alignment_input=sample.alignment_input,
    output_path=expected_path,
    number_of_shards_for_realignment=10
)
```

Getting intervals for sharding variant calling:

```python
from cpg_pipes.jobs import split_intervals
j = split_intervals.get_intervals(b=b, scatter_count=20)
```

Generate somalier pedigree fingerprints:

```python
from cpg_pipes.jobs import pedigree
fingerprint_job, fingerprint_path = pedigree.somalier_extact_job(
    b,
    sample,
    gvcf_or_cram_or_bam_path=gvcf_path,
    out_fpath=expected_path,
)
```

Check pedigree of a project:

```python
from cpg_pipes.jobs import pedigree
j, somalier_samples_path, somalier_pairs_path = pedigree.add_pedigree_jobs(
    b,
    project,
    input_path_by_sid=fingerprint_by_sid,
    web_bucket=self.pipe.web_bucket,
    web_url=self.pipe.web_url,
)
```

VQSR:

```python
from cpg_pipes.jobs import vqsr
j = vqsr.make_vqsr_jobs(
    b,
    input_vcf_or_mt_path=siteonly_vcf_path,
    gvcf_count=len(cohort.get_all_samples()),
    depends_on=inputs.get_jobs(),
    output_vcf_path=expected_path,
    use_as_annotations=True,
)
```

There other jobs available - refer to code examples.

### Pipelines

There are some full pipelines available in the `pipelines` folder, specifically:

- `pipelines/seqr_loader.py` takes CRAM/FASTQ, aligns and joint-genotypes, then annotates using Hail Query and creates an ElasticSearch index for Seqr.
- `pipelines/pedigree.py` uses Somalier and an input PED file to check pedigree and sex of CRAMs or GVCFs,
- `gatk_sv` orchestrates workflows of [GATK-SV](https://github.com/populationgenomics/gatk-sv).

### Click options

The `cpg_pipes/pipeline/cli_opts.py` module provides CLI options for Click that handle default parameters that can be used to customise a pipeline. The basic usage is the following:

```python
import click
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.pipeline import stage, Pipeline

@click.command()
@pipeline_click_options
def main(**kwargs):
    pipeline = Pipeline(
        name='my_pipeline',
        title='My pipeline',
        **kwargs,
    )
```

When calling such pipelne script from a command line, the options defined in `@pipeline_click_options` would be available and passed to the `Pipeline` initialiser:

```
  -n, --namespace [main|test|tmp] The bucket namespace to write the results to
  --analysis-project TEXT         SM project name to write the
                                  intermediate/joint-calling analysis entries
  --input-project TEXT            Only read samples that belong to the
                                  project(s). Can be set multiple times.
                                  [required]
  --first-stage TEXT              Skip previous stages and pick their expected 
                                  results if further stages depend on them
  --last-stage TEXT               Finish the pipeline after this stage
  -S, --skip-sample TEXT          Don't process specified samples. Can be set
                                  multiple times.
  -s, --only-sample TEXT          Only take these samples (can be set multiple
                                  times)
  --force-sample TEXT             Force reprocessing these samples. Can be set
                                  multiple times.
  --output-version TEXT           Suffix the outputs with this version tag.
                                  Useful for testing
  --keep-scratch / --remove-scratch
  --dry-run
  --check-smdb-seq-existence / --no-check-smdb-seq-existence
                                  Check that files in sequence.meta exist
  --skip-samples-without-first-stage-input
                                  For the first not-skipped stage, if the
                                  input for a target does notexist, just skip
                                  this target instead of failing. E.g. if the
                                  firststage is CramStage, and sequence.meta
                                  files for a sample do not exist,remove this
                                  sample instead of failing.
  --check-intermediate-existence / --no-check-intermediate-existence
                                  Within jobs, check all in-job intermediate
                                  files for possible reuse. If set to False,
                                  will overwrite all intermediates.
  --check-job-expected-outputs-existence / --no-check-job-expected-outputs-existence
                                  Before running a job, check if its input
                                  already exists. If it exists, submit a
                                  [reuse] job instead. Works nicely with
                                  --previous-batch-tsv/--previous-batch-id
                                  options.
  --update-smdb-analyses / --no-update-smdb-analyses
                                  Create analysis entries for
                                  queued/running/completed jobs
  --validate-smdb-analyses / --no-validate-smdb-analyses
                                  Validate existing analysis entries by
                                  checking if a.output exists on the bucket.
                                  Set the analysis entry to "failure" if
                                  output doesn't exist
  --previous-batch-tsv TEXT       A list of previous successful attempts from
                                  another batch, dumped from from the Batch
                                  database (the "jobs" table joined on
                                  "job_attributes"). If the intermediate
                                  output for a job exists in a previous
                                  attempt, it will be passed forward, and a
                                  [reuse] job will be submitted.
  --previous-batch-id TEXT        6-letter ID of the previous successful batch
                                  (corresponds to the directory name in the
                                  batch logs. e.g. feb0e9 in gs://cpg-seqr-
                                  main-tmp/hail/batch/feb0e9
  --source-tag TEXT               Subset found analysis to "meta={source:
                                  <source_tag>}"
  --ped-file TEXT                 PED file (will override sample-meatadata family 
                                  data if available)
  --local-dir TEXT
  --output-project TEXT           Only create ES indicies for the project(s).
                                  Can be set multiple times. Defaults to
                                  --input-projects. The name of the ES index
                                  will be suffixed with the dataset version
                                  (set by --version)
  --skip-ped-checks               Skip checking provided sex and pedigree
                                  against the inferred one
  --hc-shards-num INTEGER         Number of intervals to devide the genome for
                                  gatk HaplotypeCaller
  --use-gnarly / --no-use-gnarly  Use GnarlyGenotyper instead of GenotypeGVCFs
  --use-as-vqsr / --no-use-as-vqsr
                                  Use allele-specific annotations for VQSR
  --make-seqr-metadata / --no-make-seqr-metadata
                                  Make Seqr metadata
  --help                          Show this message and exit.
```

You can add more custom options like this:

```python
import click
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.pipeline import stage, Pipeline

@click.command()
@pipeline_click_options
@click.option('--custom-option')
def main(**kwargs):
    custom_option = kwargs.pop('custom_option')
    pipeline = Pipeline(
        name='my_pipeline',
        title='My pipeline',
        **kwargs,
    )
```

### Batch helpers.

The `cpg_pipes.hb.batch` module provides a helper function `setup_batch` to set up Hail Batch in the CPG context:

```python
from cpg_pipes.hb.batch import setup_batch
b = setup_batch('My batch')
```

It will create an instance of Batch that extends the standard Hail Batch, that records stats of added jobs and prints statistics before submission, highlighting labelled jobs, e.g.:

```
Will submit 186 jobs:
  BWA: 3 for 3 samples
  Somalier extract (CRAMs): 3 for 3 samples
  HaplotypeCaller: 3 for 3 samples
  ReblockGVCF: 3 for 3 samples
  Somalier extract (GVCFs): 3 for 3 samples
  Other jobs: 171
```

Batch instance also extend the job name if the names of a sample and a project are provided as attributes, e.g.: 

```bash
>>> j = b.new_job('My job', dict(sample='CPG196535', project='fewgenomes'))
>>> print(j.name)
fewgenomes/CPG196535: My job
```

`cpg_pipes.hb.command` provides a helper to set up a command, that can be used to add monitoring of disk space, or authentificate with Google Cloud to make gsutil work:

```python
from cpg_pipes.hb.command import wrap_command
j = b.new_job('My job')
j.command(wrap_command(
    'sleep 600',
    monitor_space=True,
    setup_gcp=True,
))
```

Will wrap the command as follows:

```
set -o pipefail
set -ex
export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

(while true; do df -h; du -sh /io; du -sh /io/batch; sleep 600; done) &

sleep 600

df -h; du -sh /io; du -sh /io/batch
```

### Reusing existing results

By default, if expected result exist, empty jobs will be submitted suffixed with ` [reuse]`. You can disable that check with `--no-check-job-expected-outputs-existence`. You can also make the pipeline trust the sample-metadata DB with `--validate-smdb-analyses` and skip checking for objects to exist on buckets.

You can also start from a specific stage with `--first-stage`, or finish on specific one with `--last-stage`.

You can use `@skipped` decorator to force skipping a stage:

```python
@skipped
@stage
class MyStage1(SampleStage):
    ...
```

`assume_results_exist=True` would also tell the code that the expected results of that stage exist, and there is no need to check bucket objects for existence:

```python
@skipped
@stage(assume_results_exist=True)
class MyStage2(SampleStage):
    ...
```

You can also force the pipeline to skip certain samples with `--skip-samples/-S`, or take only certain samples with `--only-samples/-s`, or force processing certain samples with `--force-samples`.


### Running pipelines

To start the `seqr_loader` pipeline on 2 projects `acute-care` and `perth-neuro` in the test namespace, using `seqr` as analysis project, and write index for `perth-neuro`:

```sh
python pipelines/seqr_loader.py \
-n main \
--analysis-project seqr \
--input-project acute-care \
--input-project perth-neuro \
--output-project perth-neuro \
--validate-smdb-analyses \
--check-smdb-seq-existence \
--check-intermediate-existence \
--skip-sample CPG11783 \
--skip-sample CPG13326 \
--output-version v1-0 \
--keep-scratch
```

Another useful pipeline is `pipelines/pedigree.py` to verify inferred samples 
relationship and sex against a provided PED file(s):

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
