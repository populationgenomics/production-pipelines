#!/usr/bin/env python3

"""
Test how Hail Batch attaches disks to instances. Motivation: job.storage(...) option
does not directly translate into available_storage value that you can see in job logs.
The `cpg_workflows.resources` module solves this problem by adding adjustments to the
requested storage, and this script runs some tests to determine those adjustments.
Additionally, when a larger storage is requested than available on a worker by default,
Hail Batch would attempt to attach an external disk. We want to determine at what
point it performs that attachment, as often we want to avoid extra costs incurred by
external disk, and use the default space only.
"""

from cpg_utils.config import get_config, update_dict
from cpg_utils.hail_batch import get_batch


def main():  # pylint: disable=missing-function-docstring
    update_dict(
        get_config()['workflow'],
        {
            'name': 'test_attach_disk',
            'dataset': 'fewgenomes',
            'access_level': 'test',
            'datasets': ['fewgenomes'],
        },
    )
    b = get_batch()

    def add_storage_job(
        storage: float,
        ncpu: int | None = None,
        memory: float | None = None,
    ):
        j = b.new_job(f'Disk {storage}G')
        j.storage(f'{storage}G')

        if ncpu:
            j.name += f' ncpu={ncpu}'
            j.cpu(ncpu)

        if memory:
            j.name += f' memory={memory}G'
            j.memory(f'{memory}G')

        j.command(f'sleep {60*6}')
        return j

    for _ in range(16):
        add_storage_job(storage=185 // 4, ncpu=8, memory=30.0)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
