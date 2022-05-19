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
        datasets = cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')

        output_by_target: dict[str, StageOutput | None] = dict()

        for ds_i, dataset in enumerate(datasets):
            if not dataset.get_samples():
                raise ValueError(
                    f'No active samples are found to run in the dataset {dataset}'
                )

            # Checking if all samples can be reused, queuing only one job per target:
            action_by_sid = dict()
            for sample_i, sample in enumerate(dataset.get_samples()):
                logger.info(f'{self.name}: #{sample_i + 1}/{sample}')
                action = self._get_action(sample)
                action_by_sid[sample.id] = action

            if len(set(action_by_sid.values())) == 1:
                action = list(action_by_sid.values())[0]
                if action == Action.REUSE:
                    # All stages to be reused, but adding only one reuse job
                    # (for whole dataset):
                    attrs = dataset.get_job_attrs()
                    attrs |= dict(stage=self.name, tool='[reuse]')
                    j = self.b.new_job(
                        f'{self.name} [reuse {len(dataset.get_samples())} samples]', 
                        attrs
                    )
                    inputs = self._make_inputs()
                    for _, sample in enumerate(dataset.get_samples()):
                        outputs = self.make_outputs(
                            target=sample,
                            data=self.expected_outputs(sample),
                            jobs=[j],
                        )
                        for j in outputs.jobs:
                            j.depends_on(*inputs.get_jobs(sample))
                        output_by_target[sample.target_id] = outputs
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
        output_by_target = dict()
        datasets = cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
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
