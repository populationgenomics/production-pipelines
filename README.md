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

## cpg-pipes package

The `cpg_pipes` package provides Python functions that help to design a Hail Batch powered workflow in the context of the CPG infrastructe (e.g. Google Cloud, CPG storage policies, the `analysis-runner` and the `sample-metadata` database). Specifically, the following assumptions:

* The inputs are organised in a form of projects and samples on the `sample-metadata` database. To address that, `cpg_pipes`

* `cpg_pipes.jobs` defines functions that create Hail Batch Jobs for differrent bioinformatics purposes: alignment, fastqc, deduplication, variant calling, VQSR, etc. For usage examples, see `pipelines/benchmarks/benchmark_alignment.py`, as well as other scripts in that folder.

* `cpg_pipes.hb.batch` provides a helper function `setup_batch` to set up Hail Batch in the CPG context.

Example:

```python
import click
from cpg_pipes.pipeline import \
    Pipeline, pipeline_click_options, find_stages_in_module, \
    Sample, SampleStage, StageInput, StageOutput, stage

from cpg_pipes.jobs import align

@click.command()
@pipeline_click_options
def main(**kwargs):
    # Initialize pipeline. Will automatically pass CLI options
    pipeline = Pipeline(name='my_pipeline', title='My pipeline', **kwargs)
    pipeline.set_stages(find_stages_in_module(__name__))
    pipeline.submit_batch()

@stage
class CramStage(SampleStage):
    def expected_result(self, sample: Sample):
        return f'{sample.project.get_bucket()}/cram/{sample.id}.cram'

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        expected_path = self.expected_result(self.pipe)
        job = align.bwa(b=self.pipe.b, output_path=expected_path)
        return self.make_outputs(sample, data=expected_path, jobs=[job])
```

