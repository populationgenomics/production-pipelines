"""
CPG implementation of status reporter.
"""
import os
from textwrap import dedent

from google.cloud import secretmanager
from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from ... import images, Path
from ...hb.command import wrap_command
from ...targets import Target
from ..status import (
    AnalysisStatus,
    StatusReporterError,
    StatusReporter,
)
from .smdb import SMDB, SmdbError


class CpgStatusReporter(StatusReporter):
    """
    Job status reporter. Works through creating and updating sample-metadata
    database Analysis entries.
    """

    def __init__(self, smdb: SMDB):
        self.smdb = smdb

    def add_updaters_jobs(
        self,
        b: Batch,
        output: Path | Resource | dict[str, Path | Resource],
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
    ) -> list[Job]:
        """
        Create "queued" analysis and insert "in_progress" and "completed" updater jobs.
        """
        if isinstance(output, dict):
            raise StatusReporterError(
                'SmdbStatusReporter only supports a single Path as output data.'
            )
        if isinstance(output, Resource):
            raise StatusReporterError(
                'Cannot use hail.batch.Resource objects with status reporter. '
                'Only supported single Path objects'
            )

        if not jobs:
            return []
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = self.smdb.create_analysis(
            output=str(output),
            type_=analysis_type,
            status='queued',
            sample_ids=[s.id for s in target.get_samples()],
        )
        if aid is None:
            raise StatusReporterError(
                'SmdbStatusReporter error: failed to create analysis'
            )
        # 2. Queue a job that updates the status to "in-progress"
        in_progress_j = self.add_status_updater_job(
            b,
            analysis_id=aid,
            status=AnalysisStatus.IN_PROGRESS,
            analysis_type=analysis_type,
            job_attrs=target.get_job_attrs(),
        )
        # 2. Queue a job that updates the status to "completed"
        completed_j = self.add_status_updater_job(
            b,
            analysis_id=aid,
            status=AnalysisStatus.COMPLETED,
            analysis_type=analysis_type,
            job_attrs=target.get_job_attrs(),
        )

        if prev_jobs:
            in_progress_j.depends_on(*prev_jobs)
        completed_j.depends_on(*jobs)
        return [in_progress_j, *jobs, completed_j]

    @staticmethod
    def add_status_updater_job(
        b: Batch,
        analysis_id: int,
        status: AnalysisStatus,
        analysis_type: str,
        job_attrs: dict | None = None,
    ) -> Job:
        """
        Create a Hail Batch job that updates status of analysis.
        """
        try:
            analysis_id_int = int(analysis_id)
        except ValueError:
            raise SmdbError('Analysis ID for sample-metadata must be int')

        job_name = f'Update status to {status.value}'
        if analysis_type:
            job_name += f' (for {analysis_type})'

        j = b.new_job(job_name, job_attrs)
        j.image(images.SM_IMAGE)
        cmd = dedent(
            f"""\
        cat <<EOT >> update.py
        from sample_metadata.apis import AnalysisApi
        from sample_metadata.models import AnalysisUpdateModel, AnalysisStatus
        from sample_metadata import exceptions
        import traceback
        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id={analysis_id_int},
                analysis_update_model=AnalysisUpdateModel(
                    status=AnalysisStatus('{status.value}')
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        EOT
        python update.py
        """
        )
        j.command(wrap_command(cmd, rm_leading_space=False, setup_gcp=True))
        return j
