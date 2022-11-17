"""
Metamist wrapper to report analysis progress.
"""
from enum import Enum
from textwrap import dedent
from abc import ABC, abstractmethod

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from .targets import Target
from .metamist import get_metamist, AnalysisStatus, MetamistError


class StateProvider(ABC):
    """
    Abstract pipeline state provider.
    """

    @abstractmethod
    def read_state(self, run_id) -> dict[str, dict[str, AnalysisStatus]]:
        """
        On workflow-creating time, initialise state for each stage.
        Would read state for each stage+target into a dictionary, indexed by stage ID,
        then by target ID.
        """
        pass

    @abstractmethod
    def add_status_updaters_jobs(
        self,
        b: Batch,
        outputs: dict | Path | str | None,
        stage_name: str,
        target: Target,
        analysis_type: str,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
        job_attrs: dict | None = None,
    ) -> list[Job]:
        """
        Add Hail Batch jobs that update the stage status.
        """

    @abstractmethod
    def create_status(
        self,
        outputs: dict | str | Path | None,
        status: AnalysisStatus,
        stage_name: str,
        target: Target,
        analysis_type: str,
        meta: dict | None = None,
        project_name: str = None,
    ) -> int | None:
        """
        Create status for a stage.
        """


class FSStateProvider(StateProvider):
    """
    Works through checking stage's outputs existence on the file system.
    """


class JsonFileStateProvider(StateProvider):
    """
    Works through updating a JSON file.
    """


class MetamistStateProvider(StateProvider):
    """
    Works through creating and updating Metamist's Analysis entries.
    """

    def __init__(self):
        super().__init__()

    def read_state(self):
        pass

    def add_status_updaters_jobs(
        self,
        b: Batch,
        outputs: dict | Path | str | None,
        stage_name: str,
        target: Target,
        analysis_type: str,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
        job_attrs: dict | None = None,
    ) -> list[Job]:
        """
        Create "queued" analysis and insert "in_progress" and "completed" updater jobs.
        """
        if not jobs:
            return []

        assert isinstance(outputs, str | Path | None)

        # 1. Create a "queued" analysis
        if (
            aid := self.create_status(
                outputs=outputs,
                status=AnalysisStatus.QUEUED,
                stage_name=stage_name,
                analysis_type=analysis_type,
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
            output_path=outputs if not isinstance(outputs, str | dict) else None,
        )

        if prev_jobs:
            in_progress_j.depends_on(*prev_jobs)
        completed_j.depends_on(*jobs)
        return [in_progress_j, *jobs, completed_j]

    def create_status(
        self,
        outputs: dict | str | Path | None,
        status: AnalysisStatus,
        stage_name: str,
        target: Target,
        analysis_type: str,
        meta: dict | None = None,
        dataset: str = None,
    ) -> int | None:
        """Record analysis entry"""
        assert isinstance(outputs, Path | str | None)
        return get_metamist().create_analysis(
            output=outputs,
            type_=analysis_type,
            status=status,
            sample_ids=target.get_sample_ids(),
            meta=meta,
            dataset=dataset,
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
