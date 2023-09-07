"""
Stages that implement GATK-gCNV.
"""

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_workflows.filetypes import CramPath
from cpg_workflows.jobs import gcnv
from cpg_workflows.targets import SequencingGroup, Cohort
from cpg_workflows.workflow import stage, StageInput, StageOutput
from cpg_workflows.workflow import SequencingGroupStage, CohortStage
from .. import get_batch


@stage
class PrepareIntervals(CohortStage):
    """
    Interval preparation steps that don't require the sample read counts:
    PreprocessIntervals and AnnotateIntervals.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'preprocessed': self.prefix / f'{cohort.name}.preprocessed.interval_list',
            'annotated':    self.prefix / f'{cohort.name}.annotated.tsv',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.prepare_intervals(
            get_batch(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=PrepareIntervals)
class CollectReadCounts(SequencingGroupStage):
    """
    Per-sample stage that runs CollectReadCounts to produce .counts.tsv.gz files.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'counts': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz',
            'index':  seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz.tbi',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        if seqgroup.cram is None:
            raise ValueError(f'No CRAM file found for {seqgroup}')

        jobs = gcnv.collect_read_counts(
            get_batch(),
            inputs.as_path(seqgroup.dataset, PrepareIntervals, 'preprocessed'),
            seqgroup.cram,
            self.get_job_attrs(seqgroup),
            seqgroup.dataset.prefix() / 'gcnv' / seqgroup.id,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage(required_stages=[PrepareIntervals, CollectReadCounts])
class DeterminePloidy(CohortStage):
    """
    The non-sharded cohort-wide gCNV steps after read counts have been collected:
    FilterIntervals and DetermineGermlineContigPloidy. These outputs represent
    intermediate results for the cohort as a whole, so are written to tmp_prefix.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'filtered': self.tmp_prefix / f'{cohort.name}.filtered.interval_list',
            'calls':    self.tmp_prefix / f'{cohort.name}-ploidy-calls.tar.gz',
            'model':    self.tmp_prefix / f'{cohort.name}-ploidy-model.tar.gz',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.filter_and_determine_ploidy(
            get_batch(),
            get_config()['workflow'].get('ploidy_priors'),
            inputs.as_path(cohort, PrepareIntervals, 'preprocessed'),
            inputs.as_path(cohort, PrepareIntervals, 'annotated'),
            inputs.as_path_by_target(CollectReadCounts, 'counts').values(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[PrepareIntervals, CollectReadCounts, DeterminePloidy])
class GermlineCNV(CohortStage):
    """
    The cohort-wide GermlineCNVCaller step, sharded across genome regions.
    This is separate from the DeterminePloidy stage so that the GermlineCNVCalls
    stage can pick out this stage's sharded inputs easily.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            name: self.tmp_prefix / f'{name}.tar.gz' for name in gcnv.shard_basenames()
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(cohort)

        jobs = gcnv.shard_gcnv(
            get_batch(),
            inputs.as_path(cohort, PrepareIntervals, 'annotated'),
            inputs.as_path(cohort, DeterminePloidy, 'filtered'),
            inputs.as_path(cohort, DeterminePloidy, 'calls'),
            inputs.as_path_by_target(CollectReadCounts, 'counts').values(),
            self.get_job_attrs(cohort),
            outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage(required_stages=[DeterminePloidy, GermlineCNV])
class GermlineCNVCalls(SequencingGroupStage):
    """
    Produces final individual VCF results by running PostprocessGermlineCNVCalls.
    """

    def expected_outputs(self, seqgroup: SequencingGroup) -> dict[str, Path]:
        return {
            'intervals': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.intervals.vcf.gz',
            'segments':  seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.segments.vcf.gz',
            'ratios':    seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(self, seqgroup: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        outputs = self.expected_outputs(seqgroup)

        jobs = gcnv.postprocess_calls(
            get_batch(),
            inputs.as_path(seqgroup.dataset, DeterminePloidy, 'calls'),
            inputs.as_dict(seqgroup.dataset, GermlineCNV),
            # FIXME get the sample index via sample_name.txt files instead
            seqgroup.dataset.get_sequencing_group_ids().index(seqgroup.id),
            self.get_job_attrs(seqgroup),
            outputs,
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)
