"""
Metamist wrapper to report analysis progress.
"""
import inspect
from abc import ABC, abstractmethod
from typing import Callable, Any

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch

from .targets import Target
from .metamist import get_metamist, AnalysisStatus, MetamistError


class StatusReporterError(Exception):
    """
    Error thrown by StatusReporter.
    """


class StatusReporter(ABC):
    """
    Status reporter
    """

    @abstractmethod
    def add_updaters_jobs(
        self,
        b: Batch,
        output: str,
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
    ) -> list[Job]:
        """
        Add Hail Batch jobs that update the analysis status.
        """

    @abstractmethod
    def create_analysis(
        self,
        output: str,
        analysis_type: str,
        analysis_status: str,
        target: Target,
        meta: dict | None = None,
        project_name: str = None,
    ) -> int | None:
        """
        Record analysis entry.
        """


def _calculate_size(output_path: str) -> dict[str, Any]:
    """
    Self-contained function to calculate size of an object at given path.
    @param output_path: remote path of the output file
    @return: dictionary to merge into Analysis.meta
    """
    from cloudpathlib import CloudPath

    size = CloudPath(str(output_path)).stat().st_size
    return dict(size=size)


def _update_analysis_status(
    analysis_id: int,
    new_status: AnalysisStatus,
    updater_funcs: list[Callable[[str], dict[str, Any]]] | None = None,
    output_path: str | None = None,
) -> None:
    """
    Self-contained function to update Metamist analysis entry.
    @param analysis_id: ID of Analysis entry
    @param new_status: new status to assign to the entry
    @param updater_funcs: list of functions to update the entry's metadata,
    assuming output_path as input parameter
    @param output_path: remote path of the output file, to be passed to the updaters
    """
    from sample_metadata.apis import AnalysisApi
    from sample_metadata.models import AnalysisUpdateModel
    from sample_metadata import exceptions
    from sample_metadata.model.analysis_status import AnalysisStatus as MmAnalysisStatus
    import traceback

    meta: dict[str, Any] = dict()
    if output_path and updater_funcs:
        for func in updater_funcs or []:
            meta |= func(output_path)

    aapi = AnalysisApi()
    try:
        aapi.update_analysis_status(
            analysis_id=analysis_id,
            analysis_update_model=AnalysisUpdateModel(
                status=MmAnalysisStatus(new_status.value),
                meta=meta,
            ),
        )
    except exceptions.ApiException:
        traceback.print_exc()


class MetamistStatusReporter(StatusReporter):
    """
    Job status reporter. Works through creating and updating metamist Analysis entries.
    """

    def __init__(self):
        super().__init__()

    def add_updaters_jobs(
        self,
        b: Batch,
        output: str | Path,
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
        job_attrs: dict[str, str] | None = None,
        update_analysis_meta: Callable[[str], dict] | None = None,
    ) -> list[Job]:
        """
        Create "queued" analysis and insert "in_progress" and "completed" updater jobs.
        """
        if not jobs:
            return []

        # 1. Create a "queued" analysis
        if (
            aid := self.create_analysis(
                output=str(output),
                analysis_type=analysis_type,
                analysis_status='queued',
                target=target,
                meta=meta,
            )
        ) is None:
            raise MetamistError('Failed to create analysis')

        # 2. Queue a job that updates the status to "in-progress"
        in_progress_j = self.add_status_updater_job(
            b,
            analysis_id=aid,
            status=AnalysisStatus.IN_PROGRESS,
            analysis_type=analysis_type,
            job_attrs=(job_attrs or {}) | dict(tool='metamist'),
        )
        # 2. Queue a job that updates the status to "completed"
        completed_j = self.add_status_updater_job(
            b,
            analysis_id=aid,
            status=AnalysisStatus.COMPLETED,
            analysis_type=analysis_type,
            job_attrs=(job_attrs or {}) | dict(tool='metamist'),
            output=output,
            update_analysis_meta=update_analysis_meta,
        )

        if prev_jobs:
            in_progress_j.depends_on(*prev_jobs)
        completed_j.depends_on(*jobs)
        return [in_progress_j, *jobs, completed_j]

    def create_analysis(
        self,
        output: str,
        analysis_type: str,
        analysis_status: str,
        target: Target,
        meta: dict | None = None,
        project_name: str = None,
    ) -> int | None:
        """Record analysis entry"""
        return get_metamist().create_analysis(
            output=output,
            type_=analysis_type,
            status=analysis_status,
            sample_ids=target.get_sample_ids(),
            meta=meta,
            dataset=project_name,
        )

    @staticmethod
    def add_status_updater_job(
        b: Batch,
        analysis_id: int,
        status: AnalysisStatus,
        analysis_type: str,
        job_attrs: dict | None = None,
        output: Path | str | None = None,
        update_analysis_meta: Callable[[str], dict] | None = None,
    ) -> Job:
        """
        Create a Hail Batch job that updates status of analysis. For status=COMPLETED,
        adds the size of `output` into `meta.size` if provided.
        """
        try:
            analysis_id_int = int(analysis_id)
        except ValueError:
            raise MetamistError('Analysis ID for sample-metadata must be int')

        job_name = f'Update status to {status.value}'
        if analysis_type:
            job_name += f' (for {analysis_type})'

        j = b.new_job(job_name, job_attrs)
        j.image(image_path('cpg_workflows'))

        meta_updaters_definitions = ''
        meta_updaters_funcs: list[Callable[[str], dict]] = []
        if output:
            if isinstance(output, Path):
                meta_updaters_funcs.append(_calculate_size)
            if update_analysis_meta:
                meta_updaters_funcs.append(update_analysis_meta)

            for func in meta_updaters_funcs:
                definition = inspect.getsource(func)
                if not definition.startswith('def '):
                    raise MetamistError(
                        f'Status updater must be a module-level function: {str(func)}'
                    )
                meta_updaters_definitions += definition + '\n'

        cmd = f"""
cat <<EOT >> update.py

from typing import Any, Callable
from cpg_workflows.metamist import AnalysisStatus

{meta_updaters_definitions}

{inspect.getsource(_update_analysis_status)}

{_update_analysis_status.__name__}(
    analysis_id={analysis_id_int},
    new_status=AnalysisStatus("{status.value}"),
    updater_funcs=[{', '.join(f.__name__ for f in meta_updaters_funcs)}],
    output_path={'"' + str(output) + '"' if output else 'None'},
)

EOT
python3 update.py
"""

        j.command(command(cmd, rm_leading_space=False, setup_gcp=True))
        return j
