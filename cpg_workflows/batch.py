"""
Extending the Hail's `Batch` class.
"""

import os
import tempfile
import logging
import uuid
from typing import Optional

import hailtop.batch as hb
import toml
from cpg_utils import to_path
from cpg_utils import config
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    copy_common_env,
    dataset_path,
)


_batch: Optional['Batch'] = None


def get_batch(
    name: str | None = None, default_python_image: str | None = None
) -> 'Batch':
    global _batch
    backend: hb.Backend
    if _batch is None:
        if get_config()['hail'].get('backend', 'batch') == 'local':
            logging.info('Initialising Hail Batch with local backend')
            backend = hb.LocalBackend(
                tmp_dir=tempfile.mkdtemp('batch-tmp'),
            )
        else:
            logging.info('Initialising Hail Batch with service backend')
            backend = hb.ServiceBackend(
                billing_project=get_config()['hail']['billing_project'],
                remote_tmpdir=dataset_path('batch-tmp', category='tmp'),
                token=os.environ.get('HAIL_TOKEN'),
            )
        _batch = Batch(
            name=name or get_config()['workflow'].get('name'),
            backend=backend,
            pool_label=get_config()['hail'].get('pool_label'),
            cancel_after_n_failures=get_config()['hail'].get('cancel_after_n_failures'),
            default_timeout=get_config()['hail'].get('default_timeout'),
            default_memory=get_config()['hail'].get('default_memory'),
            default_python_image=default_python_image
            or get_config()['workflow']['driver_image'],
        )
    return _batch


class Batch(hb.Batch):
    """
    Thin subclass of the Hail `Batch` class. The aim is to be able to register
    created jobs, in order to print statistics before submitting the Batch.
    """

    def __init__(
        self,
        name,
        backend,
        *args,
        pool_label=None,
        **kwargs,
    ):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry:
        self.job_by_label = dict()
        self.job_by_stage = dict()
        self.job_by_tool = dict()
        self.total_job_num = 0
        self.pool_label = pool_label
        if not get_config()['hail'].get('dry_run') and not isinstance(
            self._backend, hb.LocalBackend
        ):
            self._copy_configs_to_remote()

    def _copy_configs_to_remote(self):
        """If configs are local files, copy them to remote"""
        remote_dir = to_path(self._backend.remote_tmpdir) / 'config'
        config_path = remote_dir / (str(uuid.uuid4()) + '.toml')
        with config_path.open('w') as f:
            toml.dump(dict(get_config()), f)
        config.set_config_paths([str(config_path)])

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
        sequencing_group = attributes.get('sequencing_group')
        participant_id = attributes.get('participant_id')
        sequencing_groups: set[str] = set(attributes.get('sequencing_groups') or [])
        if sequencing_group:
            sequencing_groups.add(sequencing_group)
        part = attributes.get('part')
        label = attributes.get('label', name)
        tool = attributes.get('tool')
        if not tool and name.endswith('Dataproc cluster'):
            tool = 'hailctl dataproc'

        assert isinstance(stage, str | None)
        assert isinstance(dataset, str | None)
        assert isinstance(sequencing_group, str | None)
        assert isinstance(participant_id, str | None)
        assert isinstance(part, str | None)
        assert isinstance(label, str | None)

        name = make_job_name(
            name=name,
            sequencing_group=sequencing_group,
            participant_id=participant_id,
            dataset=dataset,
            part=part,
        )

        if label not in self.job_by_label:
            self.job_by_label[label] = {'job_n': 0, 'sequencing_groups': set()}
        self.job_by_label[label]['job_n'] += 1
        self.job_by_label[label]['sequencing_groups'] |= sequencing_groups

        if stage not in self.job_by_stage:
            self.job_by_stage[stage] = {'job_n': 0, 'sequencing_groups': set()}
        self.job_by_stage[stage]['job_n'] += 1
        self.job_by_stage[stage]['sequencing_groups'] |= sequencing_groups

        if tool not in self.job_by_tool:
            self.job_by_tool[tool] = {'job_n': 0, 'sequencing_groups': set()}
        self.job_by_tool[tool]['job_n'] += 1
        self.job_by_tool[tool]['sequencing_groups'] |= sequencing_groups

        attributes['sequencing_groups'] = list(sorted(list(sequencing_groups)))
        fixed_attrs = {k: str(v) for k, v in attributes.items()}
        return name, fixed_attrs

    def run(self, **kwargs):
        """
        Execute a batch. Overridden to print pre-submission statistics.
        """
        if not self._jobs:
            logging.error('No jobs to submit')
            return

        else:
            for job in self._jobs:
                job.name, job.attributes = self._process_attributes(
                    job.name, job.attributes
                )
                # We only have dedicated pools for preemptible machines.
                # _preemptible defaults to None, so check explicitly for False.
                if self.pool_label and job._preemptible is not False:
                    job._pool_label = self.pool_label
                copy_common_env(job)

            logging.info(f'Will submit {self.total_job_num} jobs')

            def _print_stat(_d: dict, default_label: str | None = None):
                for label, stat in _d.items():
                    label = label or default_label
                    msg = f'{stat["job_n"]} job'
                    if stat['job_n'] > 1:
                        msg += 's'
                    if (sg_count := len(stat['sequencing_groups'])) > 0:
                        msg += f' for {sg_count} sequencing group'
                        if sg_count > 1:
                            msg += 's'
                    logging.info(f'  {label}: {msg}')

            logging.info('Split by stage:')
            _print_stat(self.job_by_stage, default_label='<not in stage>')

            logging.info(f'Split by tool:')
            _print_stat(self.job_by_tool, default_label='<tool is not defined>')

        kwargs.setdefault('dry_run', get_config()['hail'].get('dry_run'))
        kwargs.setdefault(
            'delete_scratch_on_exit', get_config()['hail'].get('delete_scratch_on_exit')
        )
        if isinstance(self._backend, hb.LocalBackend):
            # Local backend does not support "wait"
            if 'wait' in kwargs:
                del kwargs['wait']
        return super().run(**kwargs)


def make_job_name(
    name: str,
    sequencing_group: str | None = None,
    participant_id: str | None = None,
    dataset: str | None = None,
    part: str | None = None,
) -> str:
    """
    Extend the descriptive job name to reflect job attributes.
    """
    if sequencing_group and participant_id:
        sequencing_group = f'{sequencing_group}/{participant_id}'
    if sequencing_group and dataset:
        name = f'{dataset}/{sequencing_group}: {name}'
    elif dataset:
        name = f'{dataset}: {name}'
    if part:
        name += f', {part}'
    return name
