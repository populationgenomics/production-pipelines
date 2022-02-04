from typing import Optional

from cpg_pipes.pipeline import Pipeline


def main():
    pipeline = Pipeline(
        analysis_project='fewgenomes',
        name='test_attach_disk',
        title='test_attach_disk',
        output_version='v0',
        namespace='test'
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
    main()
