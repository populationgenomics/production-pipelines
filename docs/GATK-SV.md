# GATK-SV

Instructions as-of 07 July 2023 for running GATK-SV via Analysis-Runner.

## Intro

[GATK-SV](https://github.com/broadinstitute/gatk-sv) is a collection of workflows for applying a number of structural
variant callers to genomic data, then refining the initial calls into a high quality joint-callset. Instead of
attempting to re-implement these workflows in Production-Pipelines, we have chosen to wrap stages of GATK-SV workflows
directly.

The [Cohort-Mode section](https://github.com/broadinstitute/gatk-sv#cohort-mode) of the README describes key GATK-SV
phases, each of which exists as a separate workflow in our implementation. This README will contain instructions for
running each of these workflows, as well as a description of required inputs.

n.b. this reflects the current state of GATK-SV adoption, but when we migrate to custom cohorts, this will be far less
complex, requiring fewer discrete steps.

## Runtime Notes

The GATK-SV workflows are run with Full permissions, as they seem to try and delete data from `tmp`, which requires
escalated permissions.

The config file `configs/gatk_sv/use_for_all_workflows.toml` contains locations of all docker images to pull from our
infrastructure, and should be included with all workflows. It also contains instruction to log analysis entries to
metamist.

## Sequencing Group Groups

To run GATK-SV we need to first identify samples which have not been processed previously. This should be done by
checking Metamist for all SGs without an Analysis entry for AnnotateVcf, the final stage of GATK-SV. All SG IDs
should be added into a TOML config file as the attribute `workflow.only_sgs`. The SGs used can span multiple projects,
so long as the relevant projects are all added to the config file as `workflow.input_datasets`, e.g.

```toml
[workflow]
input_datasets = [
    'broad-rgp',
    'acute-care',
]
only_sgs = ['CPG1', 'CPG2', 'CPG3', 'CPG4', 'CPG5']
```

There's a few considerations when building this list of SGs:

- The single sample workflow will call initial SVs on each sample individually, then analyse all samples together
  to created alike batches for processing downstream. Due to a quirk of Hail Batch ([Zulip issue](https://hail.zulipchat.com/#narrow/stream/223457-Hail-Batch-support/topic/Huge.20Network.20Setup.20Times)),
  running a massive number of SGs at once can prevent other jobs from starting. To mitigate this, we can run smaller
  batches of SGs separately, using the workflow setting `workflow.last_stages = ['GatherSampleEvidence']`, which will
  allow us to run several separate groups of SGs through the individual calling before progressing to `EvidenceQC` and
  `CreateSampleBatches`, where it is important that all SGs are run at once. If this approach is chosen, create a separate
  config TOML for each sub-batch of SGs, and run these independently.
- The amount of `tmp` space taken up by this process is wild. After completing the `GatherSampleEvidence` step for all
  samples, consider asking someone with escalated permissions in the software team to delete the data stored in
  `gs://cpg-seqr-main-tmp/cromwell/GatherSampleEvidence`, which may be tens of Terabytes of temporary files. The week-long
  standard retention policy can cost hundreds of dollars!
- The standard batch sizes for sequencing groups are set in code as 100 < 300. This can be edited in config using the
  settings `workflow.min_batch_size` and `workflow.max_batch_size`. A minimum group size of 100 is required for the
  various steps of batch processing, but The Broad implementation is sets 750 as the maximum batch size, but we have not
  experimented in our infrastructure past a max size of 120.
- One SG has repeatedly failed to run through the `GatherSampleEvidence` step, due to an elevated resource requirement.
  Samples which cause problems when generating their variant calls should be removed from the `only_sgs` variable until
  we have a way to handle these samples correctly, rather than letting it block a whole cohort.
- Some NAGIM samples still have a discrepancy between the CRAM header references and the alignment used, making them
  unusable without re-headering or realigning the files.

An example command to run the single sample workflow for a subset of SGs is:

```bash
analysis-runner \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    --config configs/genome.toml \
    --config configs/[ONLY_SGS].toml \
    --config configs/gatk_sv/use_for_all_workflows.toml \
    --dataset seqr \
    --description 'GATK-SV singlesample workflow' \
    -o gatk_sv_broad_udnaus_acute \
    --access-level full \
    python main.py gatk_sv_singlesample
```

## Separate Batch Workflows

The final stage of the `gatk_sv_singlesample` workflow is `CreateSampleBatches`. This step takes all single-SG results
and determines an optimal separation of the samples into batches for joint calling. The output of this step is a JSON
file, containing a list of all batches identified, each with:

- List of SGs to include within this batch
- Number of SGs in this batch
- Male/Female ratio of this batch
- List of Median coverage for each sample in this batch

It makes sense to evaluate the contents of this file before continuing - if the M/F ratio is significantly skewed from
1, it may make sense to add additional males or females, and repeat the single sample workflow.

If you are satisfied with the content of this file, create a new config TOML file for each of the batches identified,
following the same structure as the example above. In my proof of concept run I named these `multisam_1.toml` and
`multisam_2.toml`.

## FilterBatch

The next phase of the analysis takes the separate SG batches, and runs each separately through to FilterBatch. Each of
the separate SG groups should be run, using the config file specific to each one. The config file stop_at_filterbatch.toml
contains the instruction to stop the workflow at the correct place.

```commandline
analysis-runner \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    --config configs/genome.toml \
    --config configs/gatk_sv/multisam_1.toml \
    --config configs/gatk_sv/stop_at_filterbatch.toml \
    --config configs/gatk_sv/use_for_all_workflows.toml \
    --dataset seqr \
    --description 'GATK-SV first multisample run' \
    -o gatk_sv_broad_udnaus_acute \
    --access-level full \
    python main.py gatk_sv_multisample_1
```

## Sandwich

In the current workflow set up, the `MergeBatchSites` stage needs to be run, which will combine the `FilterBatch` data
from all separate batches into a pair of files, which can then be referenced by all batches downstream. This is done by
taking the hash of all separate batches and adding them to a new config file:

```toml
[workflow]
# this will be used when running:
# - MergeBatchSites (sandwich stage)
# - MakeCohortVCF (multisample 2)
batch_names = ['555b2843b7a3cfbdeb28c5cfbf7d2b3b1ca004_119', 'a5407cc508c8f25f960e1b52be1f31e740d289_119']
```

These hashes can be found from the output file paths generated in the previous workflow.

```commandline
analysis-runner \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    --config configs/genome.toml \
    --config configs/gatk_sv/all_sams.toml  \
    --config configs/gatk_sv/all_batch_names.toml \
    --config configs/gatk_sv/use_for_all_workflows.toml \
    --dataset seqr \
    --description 'GATK-SV first multisample run' \
    -o gatk_sv_broad_udnaus_acute \
    --access-level full \
    python main.py gatk_sv_sandwich
```

## GenotypeBatch

Once MergeBatchSites has been completed, take the paths of the two resulting output VCFs, and add to `genotypebatch.toml`.

```toml
[workflow]
only_stages = ['GenotypeBatch']
# this will be a single VCF, combined across multiple sub-cohort VCF outputs
cohort_depth_vcf = 'gs://cpg-seqr-test/gatk_sv/9beee5a0a2e112a018c97605699792d8082181_238/MergeBatchSites/cohort_depth.vcf.gz'
cohort_pesr_vcf = 'gs://cpg-seqr-main/gatk_sv/9beee5a0a2e112a018c97605699792d8082181_238/MergeBatchSites/cohort_pesr.vcf.gz'
```

Use this as a config file to finish running the `gatk_sv_multisample_1` workflow, running once for each separate batch of
SGs.

```commandline
analysis-runner \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    --config configs/genome.toml \
    --config configs/gatk_sv/multisam_1.toml \
    --config configs/gatk_sv/genotypebatch.toml \
    --config configs/gatk_sv/use_for_all_workflows.toml \
    --dataset seqr \
    --description 'GATK-SV first multisample run' \
    -o gatk_sv_broad_udnaus_acute \
    --access-level full \
    python main.py gatk_sv_multisample_1
```

## MakeCohortVCF + AnnotateVCF

The final stage of GATK-SV takes all the jointly genotyped data from the separate batches and combines them to create
a single VCF file, which is annotated with functional and AF annotations. This runs on all samples together, using the
`only_sgs` variable from the Single sample workflow.

```commandline
analysis-runner \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    --config configs/genome.toml \
    --config configs/gatk_sv/all_sams.toml \
    --config configs/gatk_sv/all_batch_names.toml \
    --config configs/gatk_sv/use_for_all_workflows.toml \
    --dataset seqr \
    --description 'GATK-SV final multisample run' \
    -o gatk_sv_broad_udnaus_acute \
    --access-level full \
    python main.py gatk_sv_multisample_2
```

## Even More Batches

In theory, every sub-batch ID previously analysed can be combined at the Sandwich stage, and for the gatk_sv_multisample_2
workflow, to create a union joint-call across every SG analysed so far. This might become increasingly expensive to
re-compute.
