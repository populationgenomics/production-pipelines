"""
CPG implementation of status reporter.
"""

import os
from textwrap import dedent

from cpg_utils.hail_batch import image_path
from google.cloud import secretmanager
from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from ... import Path
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
    database Analysis entries. It is also able to send notifications to Slack.
    To enable that, create a channel, set `CPG_SLACK_CHANNEL` and `CPG_SLACK_TOKEN` 
    environment  variables, and add "Seqr Loader" app into a channel with:

    /invite @Seqr Loader
    """

    def __init__(
        self, 
        smdb: SMDB, 
        slack_channel: str | None = None,
    ):
        super().__init__()
        self.smdb = smdb
        self.slack_channel = slack_channel or os.environ.get('CPG_SLACK_CHANNEL')
        self.slack_token = os.environ.get('CPG_SLACK_TOKEN')
        if self.slack_channel and not self.slack_token:
            project_id = 'seqr-308602'
            secret_name = 'slack-seqr-loader-token'
            slack_token_secret = (
                f'projects/{project_id}/secrets/{secret_name}/versions/latest'
            )
            secret_manager = secretmanager.SecretManagerServiceClient()
            # noinspection PyTypeChecker
            response = secret_manager.access_secret_version(
                request={'name': slack_token_secret}
            )
            self.slack_token = response.payload.data.decode('UTF-8')

    def add_updaters_jobs(
        self,
        b: Batch,
        output: str | Path | Resource | dict[str, Path | Resource],
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        prev_jobs: list[Job] | None = None,
        meta: dict | None = None,
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
        if (aid := self.create_analysis(
            output=str(output),
            analysis_type=analysis_type,
            analysis_status='queued',
            target=target,
            meta=meta,
        )) is None:
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
            output=output if isinstance(output, Path) else None,
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
    ) -> int | None:
        """Record analysis entry"""
        return self.smdb.create_analysis(
            output=output,
            type_=analysis_type,
            status=analysis_status,
            sample_ids=target.get_sample_ids(),
            meta=meta,
        )

    @staticmethod
    def add_status_updater_job(
        b: Batch,
        analysis_id: int,
        status: AnalysisStatus,
        analysis_type: str,
        job_attrs: dict | None = None,
        output: Path | None = None,
    ) -> Job:
        """
        Create a Hail Batch job that updates status of analysis. For status=COMPLETED,
        adds the size of `output` into `meta.size` if provided.
        """
        try:
            analysis_id_int = int(analysis_id)
        except ValueError:
            raise SmdbError('Analysis ID for sample-metadata must be int')

        job_name = f'Update status to {status.value}'
        if analysis_type:
            job_name += f' (for {analysis_type})'

        j = b.new_job(job_name, job_attrs)
        j.image(image_path('sm-api'))

        calc_size_cmd = None
        if output:
            calc_size_cmd = f"""
        from cloudpathlib import CloudPath
        meta['size'] = CloudPath('{str(output)}').stat().st_size
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
        python update.py
        """
        )
        j.command(wrap_command(cmd, rm_leading_space=False, setup_gcp=True))
        return j
