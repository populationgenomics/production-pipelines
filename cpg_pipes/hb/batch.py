"""
Extending the Hail's `Batch` class.
"""

import logging
import os

import hailtop.batch as hb
from cloudpathlib import CloudPath
from cpg_utils import to_path
from cpg_utils import config
from cpg_utils.config import get_config
from cpg_utils.hail_batch import copy_common_env, remote_tmpdir


logger = logging.getLogger(__file__)


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
        attributes: dict | None = None,
    ) -> tuple[str, dict[str, str]]:
        """
        Use job attributes to make the job name more descriptive, and add
        labels for Batch pre-submission stats.
        """
        if not name:
            raise ValueError('Error: job name must be defined')

        self.total_job_num += 1

        attributes = attributes or {}
        stage = attributes.get('stage')
        dataset = attributes.get('dataset')
        sample = attributes.get('sample')
        participant_id = attributes.get('participant_id')
        samples: set[str] = set(attributes.get('samples') or [])
        if sample:
            samples.add(sample)
        part = attributes.get('part')
        label = attributes.get('label', name)
        tool = attributes.get('tool')
        reuse = attributes.get('reuse', False)
        if reuse and not tool:
            tool = '[reuse]'

        assert isinstance(stage, str | None)
        assert isinstance(dataset, str | None)
        assert isinstance(sample, str | None)
        assert isinstance(participant_id, str | None)
        assert isinstance(part, str | None)
        assert isinstance(label, str | None)
        assert isinstance(reuse, bool)

        name = make_job_name(
            name=name,
            sample=sample,
            participant_id=participant_id,
            dataset=dataset,
            part=part,
            reuse=reuse,
        )

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

    def run(self, **kwargs):
        """
        Execute a batch. Overridden to print pre-submission statistics.
        """
        if not self._jobs:
            logger.error('No jobs to submit')
            return

        else:
            for job in self._jobs:
                job.name, job.attributes = self._process_attributes(
                    job.name, job.attributes
                )
                if self.pool_label:
                    job._pool_label = self.pool_label
                copy_common_env(job)

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

            logger.info(f'Split by tool:')
            _print_stat(self.job_by_tool, default_label='<tool is not defined>')

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
    access_level = get_config()['workflow']['access_level']
    pool_label = get_config()['hail'].get('pool_label')

    bucket = remote_tmpdir(f'cpg-{dataset}-hail')
    if access_level == 'test':
        bucket = remote_tmpdir(f'cpg-{dataset}-hail/test')

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
    participant_id: str | None = None,
    dataset: str | None = None,
    part: str | None = None,
    reuse: bool = False,
) -> str:
    """
    Extend the descriptive job name to reflect job attributes.
    """
    if sample and participant_id:
        sample = f'{sample}/{participant_id}'
    if sample and dataset:
        name = f'{dataset}/{sample}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    if part:
        name += f', {part}'
    if reuse:
        name += ' [reuse]'
    return name
