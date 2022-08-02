"""
Extending the Hail's `Batch` class.
"""

import logging
import os
from typing import TypedDict

import hailtop.batch as hb
from cloudpathlib import CloudPath
from cpg_utils import to_path
from cpg_utils import config
from cpg_utils.config import get_config
from cpg_utils.hail_batch import copy_common_env, remote_tmpdir
from hailtop.batch.job import PythonJob, BashJob


logger = logging.getLogger(__file__)


class JobAttributes(TypedDict, total=False):
    """
    Job attributes specification.
    """

    sample: str
    dataset: str
    samples: list[str]
    datasets: list[str]
    part: str
    label: str
    stage: str
    tool: str
    reuse: bool


class RegisteringBatch(hb.Batch):
    """
    Thin subclass of the Hail `Batch` class. The aim is to be able to register
    created jobs, in order to print statistics before submitting the Batch.
    """

    def __init__(self, name, backend, *args, pool_label=None, **kwargs):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry:
        self.job_by_label = dict()
        self.job_by_stage = dict()
        self.job_by_tool = dict()
        self.total_job_num = 0
        self.pool_label = pool_label
        self._copy_configs_to_remote()

    def _copy_configs_to_remote(self):
        """If configs are local files, copy them to remote"""
        remote_dir = to_path(self._backend.remote_tmpdir) / 'config'
        remote_paths = []
        # noinspection PyProtectedMember
        for path in config._config_paths:
            path = to_path(path)
            if isinstance(path, CloudPath):
                remote_paths.append(str(path))
            else:
                remote_path = remote_dir / path.name
                with path.open() as inp, remote_path.open('w') as out:
                    out.write(inp.read())
                remote_paths.append(str(remote_path))
        config.set_config_paths(remote_paths)
        os.environ['CPG_CONFIG_PATH'] = ','.join(remote_paths)

    def _process_attributes(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
    ) -> tuple[str, dict[str, str]]:
        """
        Use job attributes to make the job name more descriptive, and add
        labels for Batch pre-submission stats.
        """
        if not name:
            raise ValueError('Error: job name must be defined')

        self.total_job_num += 1

        attributes = attributes or JobAttributes()
        stage = attributes.get('stage')
        dataset = attributes.get('dataset')
        sample = attributes.get('sample')
        samples: set[str] = set(attributes.get('samples') or [])
        if sample:
            samples.add(sample)
        part = attributes.get('part')
        label = attributes.get('label', name)
        tool = attributes.get('tool')
        reuse = attributes.get('reuse', False)
        if reuse and not tool:
            tool = '[reuse]'

        name = make_job_name(name, sample, dataset, part, reuse)

        if label not in self.job_by_label:
            self.job_by_label[label] = {'job_n': 0, 'samples': set()}
        self.job_by_label[label]['job_n'] += 1
        self.job_by_label[label]['samples'] |= samples

        if stage not in self.job_by_stage:
            self.job_by_stage[stage] = {'job_n': 0, 'samples': set()}
        self.job_by_stage[stage]['job_n'] += 1
        self.job_by_stage[stage]['samples'] |= samples

        if tool not in self.job_by_tool:
            self.job_by_tool[tool] = {'job_n': 0, 'samples': set()}
        self.job_by_tool[tool]['job_n'] += 1
        self.job_by_tool[tool]['samples'] |= samples

        attributes['samples'] = list(sorted(list(samples)))
        fixed_attrs = {k: str(v) for k, v in attributes.items()}
        return name, fixed_attrs

    def new_job(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
        **kwargs,
    ) -> BashJob:
        """
        Wrapper around `new_job()` that processes job attributes.
        """
        name, fixed_attributes = self._process_attributes(name, attributes)
        j = super().new_job(name, attributes=fixed_attributes, **kwargs)
        if self.pool_label:
            j._pool_label = self.pool_label
        copy_common_env(j)
        return j

    def new_python_job(
        self,
        name: str | None = None,
        attributes: JobAttributes | None = None,
    ) -> PythonJob:
        """
        Wrapper around `new_python_job()` that processes job attributes.
        """
        name, fixed_attributes = self._process_attributes(name, attributes)
        j = super().new_python_job(name, attributes=fixed_attributes)
        if self.pool_label:
            j._pool_label = self.pool_label
        copy_common_env(j)
        return j

    def run(self, **kwargs):
        """
        Execute a batch. Overridden to print pre-submission statistics.
        """
        logger.info(f'Will submit {self.total_job_num} jobs')

        def _print_stat(_d: dict, default_label: str | None = None):
            for label, stat in _d.items():
                label = label or default_label
                msg = f'{stat["job_n"]} job'
                if stat['job_n'] > 1:
                    msg += 's'
                if len(stat['samples']) > 0:
                    msg += f' for {len(stat["samples"])} sample'
                    if len(stat['samples']) > 1:
                        msg += 's'
                logger.info(f'  {label}: {msg}')

        logger.info('Split by stage:')
        _print_stat(self.job_by_stage, default_label='<not in stage>')

        logger.info('Split by label:')
        _print_stat(self.job_by_stage, default_label='<no label>')

        logger.info(f'Split by tool:')
        _print_stat(self.job_by_tool, default_label='<tool is not defined>')

        if get_config().get('dry_run', False):
            return

        return super().run(
            dry_run=get_config()['hail'].get('dry_run', False),
            delete_scratch_on_exit=not get_config()['hail'].get('keep_scratch', True),
            **kwargs,
        )


def setup_batch(description: str) -> RegisteringBatch:
    """
    Wrapper around the initialisation of a Hail Batch object.

    @param description: descriptive name of the Batch (will be displayed in the GUI)
    """
    billing_project = get_config()['hail']['billing_project']
    dataset = get_config()['workflow']['dataset']
    pool_label = get_config()['hail'].get('pool_label')
    bucket = remote_tmpdir(f'cpg-{dataset}-hail')

    logger.info(
        f'Starting Hail Batch with the project {billing_project}'
        f', bucket {bucket}' + (f', pool label {pool_label}' if pool_label else '')
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=bucket,
        token=os.environ.get('HAIL_TOKEN'),
    )
    return RegisteringBatch(
        name=description,
        backend=backend,
        pool_label=pool_label,
    )


def make_job_name(
    name: str,
    sample: str | None = None,
    dataset: str | None = None,
    part: str | None = None,
    reuse: bool = False,
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
    if reuse:
        name += ' [reuse]'
    return name
