"""
Metamist wrapper to report analysis progress.
"""

from textwrap import dedent
from abc import ABC, abstractmethod

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from .targets import Target
from .metamist import get_metamist, AnalysisStatus, MetamistError


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


class MetamistStatusReporter(StatusReporter):
    """
    Job status reporter. Works through creating and updating metamist Analysis entries.
    """

    def __init__(self):
        super().__init__()

    def add_updaters_jobs(
        self,
        b: Batch,
        output: str,
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
        job_attrs: dict[str, str] | None = None,
    ) -> list[Job]:
        """
        Create "queued" analysis and insert "in_progress" and "completed" updater jobs.
        """
        if not jobs:
            return []

        # 1. Create a "queued" analysis
        if (
            aid := self.create_analysis(
                output=output,
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
            output_path=output if not isinstance(output, str | dict) else None,
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
        output_path: Path | None = None,
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

        calc_size_cmd = None
        if output_path:
            calc_size_cmd = f"""
        from cloudpathlib import CloudPath
        meta['size'] = CloudPath('{str(output_path)}').stat().st_size
        """
        cmd = dedent(
            f"""\
        cat <<EOT >> update.py
        from sample_metadata.apis import AnalysisApi
        from sample_metadata.models import AnalysisUpdateModel, AnalysisStatus
        from sample_metadata import exceptions
        import traceback

        meta = dict()
        {calc_size_cmd}

        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id={analysis_id_int},
                analysis_update_model=AnalysisUpdateModel(
                    status=AnalysisStatus('{status.value}'),
                    meta=meta,
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        EOT
        python3 update.py
        """
        )
        j.command(command(cmd, rm_leading_space=False, setup_gcp=True))
        return j
