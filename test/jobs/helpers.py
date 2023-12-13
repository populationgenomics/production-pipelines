from hailtop.batch.job import Job
from hailtop.batch.resource import ResourceFile


def get_path_from_resource_file(file: ResourceFile) -> str:
    path = str(getattr(file, '_input_path', ''))
    if not path:
        raise ValueError(f"ResourceFile '{file}' must have an '_input_path' attribute")
    return path


def get_command_str(job: Job) -> str:
    return ''.join(getattr(job, '_command', []))
