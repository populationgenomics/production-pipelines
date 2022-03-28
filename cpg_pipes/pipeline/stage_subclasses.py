"""
Subclasses of `Stage` specific for a target subclass.
"""
import logging
from abc import ABC, abstractmethod

from .pipeline import Pipeline, Stage, ExpectedResultT, StageInput, StageOutput
from .targets import Sample, Cohort, Dataset

logger = logging.getLogger(__file__)


class SampleStage(Stage[Sample], ABC):
    """
    Sample-level stage.
    """

    @abstractmethod
    def expected_outputs(self, sample: Sample) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Pplug in stage into the pipeline.
        """
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for ds_i, ds in enumerate(datasets):
            logger.info(f'{self.name}: #{ds_i} {ds}')
            if not ds.get_samples():
                raise ValueError(
                    f'No active samples are found to run in the dataset {ds.name}'
                )
            for sample_i, sample in enumerate(ds.get_samples()):
                logger.info(f'{self.name}: #{sample_i}/{sample}')
                sample_result = self._queue_jobs_with_checks(sample)
                output_by_target[sample.target_id] = sample_result
            logger.info('-#-#-#-')
        return output_by_target


class DatasetStage(Stage, ABC):
    """
    Dataset-level stage
    """

    @abstractmethod
    def expected_outputs(self, dataset: Dataset) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Pplug in stage into the pipeline.
        """
        output_by_target = dict()
        datasets = pipeline.cohort.get_datasets()
        if not datasets:
            raise ValueError('No active datasets are found to run')
        for _, ds in enumerate(datasets):
            output_by_target[ds.target_id] = self._queue_jobs_with_checks(ds)
        return output_by_target


class CohortStage(Stage, ABC):
    """
    Entire cohort level stage
    """

    @abstractmethod
    def expected_outputs(self, cohort: Cohort) -> ExpectedResultT:
        """
        to_path(s) to files that the stage is epxected to generate for the `target`.
        Used within the stage to pass the output paths to commands, as well as
        by the pipeline to get expected paths when the stage is skipped and
        didn't return a `StageDeps` object from `add_jobs()`.
        """

    @abstractmethod
    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        """
        Override to add Hail Batch jobs.
        """
        pass

    def add_to_the_pipeline(self, pipeline: Pipeline) -> dict[str, StageOutput]:
        """
        Override to plug in stage into the pipeline.
        """
        return {
            pipeline.cohort.target_id:
                self._queue_jobs_with_checks(pipeline.cohort)
        }
