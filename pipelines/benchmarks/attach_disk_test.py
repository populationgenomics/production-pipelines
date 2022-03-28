"""
Test how Hail Batch attaches disks to instances
"""

from typing import Optional

from cpg_pipes import Namespace
from cpg_pipes.pipeline import create_pipeline


def main():  # pylint: disable=missing-function-docstring
    pipeline = create_pipeline(
        analysis_dataset='fewgenomes',
        name='test_attach_disk',
        description='test_attach_disk',
        version='v0',
        namespace=Namespace.TEST,
    )
    b = pipeline.b

    def add_storage_job(
        storage: float, 
        ncpu: Optional[int] = None,
        memory: Optional[float] = None,
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

    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
