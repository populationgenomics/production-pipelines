"""
Extending the Hail's `Batch` class.
"""

from hail.backend.service_backend import ServiceBackend
from hail.utils.java import Env


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


_override_revision = None


class OverrideServiceBackend(ServiceBackend):
    @property
    def jar_spec(self) -> dict:
        return {'type': 'git_revision', 'value': _override_revision}


def override_jar_spec(revision: str) -> None:
    global _override_revision
    _override_revision = revision
    backend = Env.backend()
    if isinstance(backend, ServiceBackend):
        backend.__class__ = OverrideServiceBackend
