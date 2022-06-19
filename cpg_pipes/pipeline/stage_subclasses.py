"""
Subclasses of Stage class specific for a Target subclass.
"""

import logging
from abc import ABC, abstractmethod

from .pipeline import Stage, ExpectedResultT, StageInput, StageOutput, Action
from ..targets import Sample, Cohort, Dataset

logger = logging.getLogger(__file__)


class SampleStage(Stage[Sample], ABC):
    """
    Sample-level stage.
    """

    @abstractmethod
    def expected_outputs(self, sample: Sample) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug in stage into the pipeline.
        """
        output_by_target: dict[str, StageOutput | None] = dict()

        if not (datasets := cohort.get_datasets()):
            logger.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.datasets` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`'
            )
            return output_by_target
        if not cohort.get_samples():
            logger.warning(
                f'{len(cohort.get_samples())}/'
                f'{len(cohort.get_samples(only_active=False))} '
                f'usable (active=True) samples found. Check logs above for '
                f'possible reasons samples were skipped (e.g. all samples ignored '
                f'via `workflow.skip_samples` in config, or they all missing stage '
                f'inputs and `workflow.skip_samples_with_missing_input=true` is set)'
            )
            return output_by_target

        for ds_i, dataset in enumerate(datasets):
            if not dataset.get_samples():
                logger.warning(
                    f'{dataset}: '
                    f'{len(dataset.get_samples())}/'
                    f'{len(dataset.get_samples(only_active=False))} '
                    f'usable (active=True) samples found. Check logs above for '
                    f'possible reasons samples were skipped (e.g. all samples ignored '
                    f'via `workflow.skip_samples` in config, or they all missing stage '
                    f'inputs and `workflow.skip_samples_with_missing_input=true` is set)'
                )
                continue

            # Checking if all samples can be reused, queuing only one job per target:
            action_by_sid = dict()
            for sample_i, sample in enumerate(dataset.get_samples()):
                logger.info(f'{self.name}: #{sample_i + 1}/{sample}')
                action = self._get_action(sample)
                action_by_sid[sample.id] = action
                if action == Action.REUSE:
                    if self.analysis_type and self.status_reporter:
                        self.status_reporter.create_analysis(
                            output=str(self.expected_outputs(sample)),
                            analysis_type=self.analysis_type,
                            analysis_status='completed',
                            target=sample,
                            meta=sample.get_job_attrs(),
                        )

            if len(set(action_by_sid.values())) == 1:
                action = list(action_by_sid.values())[0]
                if action == Action.REUSE:
                    for _, sample in enumerate(dataset.get_samples()):
                        output_by_target[sample.target_id] = self.make_outputs(
                            target=sample,
                            data=self.expected_outputs(sample),
                        )
                    continue

            # Some samples can't be reused, queuing each sample:
            logger.info(f'{self.name}: #{ds_i + 1} {dataset}')
            for sample_i, sample in enumerate(dataset.get_samples()):
                logger.info(f'{self.name}: #{sample_i + 1}/{sample}')
                output_by_target[sample.target_id] = self._queue_jobs_with_checks(
                    sample, action=action_by_sid[sample.id]
                )
            logger.info('-#-#-#-')

        return output_by_target


class DatasetStage(Stage, ABC):
    """
    Dataset-level stage
    """

    @abstractmethod
    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Plug in stage into the pipeline.
        """
        output_by_target: dict[str, StageOutput | None] = dict()
        if not (datasets := cohort.get_datasets()):
            logger.warning(
                f'{len(cohort.get_datasets())}/'
                f'{len(cohort.get_datasets(only_active=False))} '
                f'usable (active=True) datasets found in the cohort. Check that '
                f'`workflow.datasets` is provided, and not all datasets are skipped '
                f'via workflow.skip_datasets`'
            )
            return output_by_target
        for _, dataset in enumerate(datasets):
            output_by_target[dataset.target_id] = self._queue_jobs_with_checks(dataset)
        return output_by_target


class CohortStage(Stage, ABC):
    """
    Cohort-level stage (all datasets of a pipeline run).
    """

    @abstractmethod
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        """
        Override to declare expected output paths.
        """

    @abstractmethod
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def queue_for_cohort(self, cohort: Cohort) -> dict[str, StageOutput | None]:
        """
        Override to plug in stage into the pipeline.
        """
        return {cohort.target_id: self._queue_jobs_with_checks(cohort)}
