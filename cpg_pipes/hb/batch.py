"""
Extending the Hail's `Batch` class.
"""

import logging
import os
from typing import TypedDict

import hailtop.batch as hb
from hailtop.batch.job import Job, PythonJob, BashJob

from .. import Path

logger = logging.getLogger(__file__)


class JobAttributes(TypedDict, total=False):
    """
    Job attributes specification.
    """

    sample: str
    dataset: str
    samples: str
    part: str
    label: str
    stage: str
    tool: str


class RegisteringBatch(hb.Batch):
    """
    Thin subclass of the Hail `Batch` class. The aim is to be able to register
    created jobs, in order to print statistics before submitting the Batch.
    """

    def __init__(
        self, name, backend, *args, use_private_pool=False, **kwargs
    ):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry:
        self.job_by_label = dict()
        self.job_by_stage = dict()
        self.job_by_tool = dict()
        self.total_job_num = 0
        self.use_private_pool = use_private_pool

    def _process_attributes(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
    ) -> str:
        """
        Use job attributes to make the job name more descriptibe, and add
        labels for Batch pre-submission stats.
        """
        if not name:
            raise ValueError('Error: job name must be defined')

        self.total_job_num += 1

        attributes = attributes or JobAttributes()
        stage = attributes.get('stage')
        dataset = attributes.get('dataset')
        sample = attributes.get('sample')
        sample_list = attributes.get('samples')
        part = attributes.get('part')
        label = attributes.get('label', name)
        tool = attributes.get('tool')

        name = make_job_name(name, sample, dataset, part)

        samples = []
        if sample_list:
            samples = sample_list.split(',')
        elif sample:
            samples = [sample]

        if label not in self.job_by_label:
            self.job_by_label[label] = {'job_n': 0, 'samples': set()}
        self.job_by_label[label]['job_n'] += 1
        self.job_by_label[label]['samples'] |= set(samples)

        if stage not in self.job_by_stage:
            self.job_by_stage[stage] = {'job_n': 0, 'samples': set()}
        self.job_by_stage[stage]['job_n'] += 1
        self.job_by_stage[stage]['samples'] |= set(samples)
        
        if tool not in self.job_by_tool:
            self.job_by_tool[tool] = {'job_n': 0, 'samples': set()}
        self.job_by_tool[tool]['job_n'] += 1
        self.job_by_tool[tool]['samples'] |= set(samples)
        
        return name

    def new_python_job(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
    ) -> PythonJob:
        """
        Wrapper around `new_python_job()` that processes job attributes.
        """
        name = self._process_attributes(name, attributes)
        return super().new_python_job(name, attributes=attributes)

    def new_job(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
        **kwargs,
    ) -> BashJob:
        """
        Wrapper around `new_job()` that processes job attributes.
        """
        name = self._process_attributes(name, attributes)
        j = super().new_job(name, attributes=attributes)
        if self.use_private_pool:
            j.attributes['use_private_pool'] = 'TRUE'
        return j

    def run(self, **kwargs):
        """
        Execute a batch. Overridden to print pre-submission statistics.
        """
        logger.info(f'Will submit {self.total_job_num} jobs. Split by stage:')
        for stage, stat in self.job_by_stage.items():
            logger.info(
                f'  {stage or "<not in stage>"}: '
                f'{stat["job_n"]} jobs for {len(stat["samples"])} samples'
            )

        logger.info(f'Split by label:')
        for label, stat in self.job_by_label.items():
            logger.info(
                f'  {label}: {stat["job_n"]} jobs for {len(stat["samples"])} samples'
            )

        logger.info(f'Split by tool:')
        for tool, stat in self.job_by_tool.items():
            logger.info(
                f'  {tool or "<tool is not defined>"}: '
                f'{stat["job_n"]} jobs for {len(stat["samples"])} samples'
            )

        return super().run(**kwargs)


def setup_batch(
    description: str,
    billing_project: str | None = None,
    hail_bucket: Path | None = None,
    use_private_pool: bool = False,
) -> RegisteringBatch:
    """
    Wrapper around the initialization of a Hail Batch object.
    Handles setting the temporary bucket and the billing project.

    @param description: descriptive name of the Batch (will be displayed in the GUI)
    @param billing_project: Hail billing project name
    @param hail_bucket: bucket for Hail Batch intermediate files.
    @param use_private_pool: submit jobs to the private pool
    """
    billing_project = get_billing_project(billing_project)
    hail_bucket = get_hail_bucket(hail_bucket)

    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=hail_bucket,
        token=os.environ.get('HAIL_TOKEN'),
    )
    return RegisteringBatch(
        name=description, 
        backend=backend, 
        use_private_pool=use_private_pool
    )


def get_hail_bucket(hail_bucket: str | Path | None = None) -> str:
    """
    Get Hail bucket.
    """
    if not hail_bucket:
        hail_bucket = os.getenv('HAIL_BUCKET')
        if not hail_bucket:
            raise ValueError(
                'Either the hail_bucket parameter, or the HAIL_BUCKET '
                'environment variable must be set'
            )
    return str(hail_bucket)


def get_billing_project(billing_project: str | None = None) -> str:
    """
    Get Hail billing project.
    """

    billing_project = billing_project or os.getenv('HAIL_BILLING_PROJECT')
    if not billing_project:
        raise ValueError(
            'Either the billing_project parameter, or the HAIL_BILLING_PROJECT '
            'environment variable must be set'
        )
    return billing_project


def make_job_name(
    name: str,
    sample: str | None = None,
    dataset: str | None = None,
    part: str | None = None,
) -> str:
    """
    Extend the descriptive job name to reflect job attributes.
    """
    if sample and dataset:
        name = f'{dataset}/{sample}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    if part:
        name += f', {part}'
    return name


def hail_query_env(
    j: Job, hail_billing_project: str, hail_bucket: str | Path | None = None
):
    """
    Setup environment to run Hail Query Service backend script.
    """
    j.env('HAIL_BILLING_PROJECT', hail_billing_project)
    j.env('HAIL_BUCKET', get_hail_bucket(hail_bucket))
