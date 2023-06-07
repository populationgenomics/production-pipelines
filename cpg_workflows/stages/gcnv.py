"""
Stages that run gCNV.
"""

from cpg_utils.config import get_config
from cpg_workflows.jobs import gatk
from cpg_workflows.workflow import Sample, stage, StageInput, StageOutput, Dataset, DatasetStage, SampleStage
from cpg_workflows.stages.align import Align
from .. import get_batch


@stage(required_stages=Align)
class PrepareIntervals(DatasetStage):
    def _intervals_file(self, dataset):
        h = dataset.alignment_inputs_hash()
        return dataset.prefix() / 'gcnv' / h / f'{dataset.name}.interval_list'

    def expected_outputs(self, dataset):
        return {'intervals': self._intervals_file(dataset)}

    def queue_jobs(self, dataset, inputs):
        jobs = gatk.preprocess_intervals(
            get_batch(),
            get_config()['workflow'].get('intervals_path'),
            self.get_job_attrs(dataset),
            self._intervals_file(dataset)
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)


@stage(required_stages=PrepareIntervals)
class CollectReadCounts(SampleStage):
    """
    Per-sample stage that produces .counts.hdf5 files.
    """

    def _counts_file(self, sample):
        return sample.dataset.prefix() / 'gcnv' / f'{sample.id}.counts.hdf5'

    def expected_outputs(self, sample):
        return {'gcnv_counts': self._counts_file(sample)}

    def queue_jobs(self, sample, inputs):
        job_attrs = self.get_job_attrs(sample)
        jobs = gatk.collect_read_counts(
            get_batch(),
            sample,
            inputs.as_path(sample.dataset, PrepareIntervals, 'intervals'),
            job_attrs,
            self._counts_file(sample)
        )
        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)


@stage(required_stages=CollectReadCounts)
class DeterminePloidy(DatasetStage):
    def _file_dir(self, dataset):
        h = dataset.alignment_inputs_hash()
        return dataset.prefix() / 'gcnv' / h

    def expected_outputs(self, dataset):
        return {
            'model': self._file_dir(dataset) / f'{dataset.name}-contig-ploidy-model.tar.gz',
            'calls': self._file_dir(dataset) / f'{dataset.name}-contig-ploidy-calls.tar.gz'
        }

    def queue_jobs(self, dataset, inputs):
        jobs = gatk.determine_ploidy(
            get_batch(),
            dataset.name,
            get_config()['workflow'].get('ploidy_priors'),
            [inputs.as_path(s, CollectReadCounts, 'gcnv_counts') for s in dataset.get_samples()],
            self.get_job_attrs(dataset),
            self._file_dir(dataset) / f'{dataset.name}-contig-ploidy'
        )
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
