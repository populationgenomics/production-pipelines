"""
Record workflow progress.
"""
import logging
from textwrap import dedent
from abc import ABC, abstractmethod

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch, Resource

from . import get_cohort
from .targets import Target, Cohort
from .metamist import get_metamist, AnalysisStatus, MetamistError, AnalysisType


class StateProvider(ABC):
    """
    Abstract pipeline state provider.
    """

    @abstractmethod
    def read_state(
        self, cohort: Cohort, run_id: str
    ) -> dict[str, dict[str, AnalysisStatus]]:
        """
        On workflow-creating time, initialise state for each stage.
        Would read state for each stage+target into a dictionary, indexed by stage ID,
        then by target ID.
        """
        pass

    @abstractmethod
    def record_status(
        self,
        outputs: dict | str | Path | None,
        status: AnalysisStatus,
        stage_name: str,
        target: Target,
        analysis_type: str,
        meta: dict | None = None,
        dataset: str | None = None,
        main_output_key: str | None = None,
    ) -> int:
        """
        Record status of a stage
        """

    def wrap_jobs_with_status_updaters(
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
        main_output_key: str | None = None,
    ) -> list[Job]:
        """
        Record QUEUED status for a stage, and insert jobs that update status to
        IN_PROGRESS and COMPLETED.
        """
        if not jobs:
            return []

        entry_id: int = self.record_status(
            outputs=outputs,
            status=AnalysisStatus.QUEUED,
            stage_name=stage_name,
            analysis_type=analysis_type,
            target=target,
            meta=meta,
            main_output_key=main_output_key,
        )

        # 2. Queue a job that updates the status to IN_PROGRESS
        in_progress_j = self.add_status_updater_job(
            b,
            entry_id=entry_id,
            status=AnalysisStatus.IN_PROGRESS,
            job_attrs=job_attrs,
        )
        # 2. Queue a job that updates the status to COMPLETED
        completed_j = self.add_status_updater_job(
            b,
            entry_id=entry_id,
            status=AnalysisStatus.COMPLETED,
            job_attrs=job_attrs,
        )

        in_progress_j.depends_on(*(prev_jobs or []))
        completed_j.depends_on(*jobs)
        return [in_progress_j, *jobs, completed_j]

    def add_status_updater_job(
        self,
        b: Batch,
        entry_id: int,
        status: AnalysisStatus,
        analysis_type: str | None = None,
        job_attrs: dict | None = None,
        outputs: dict | str | Path | None = None,
        main_output_key: str | None = None,
    ) -> Job:
        """
        Create a Hail Batch job that updates status of analysis. For status=COMPLETED,
        adds the size of `output` into `meta.size` if provided.
        """
        job_name = f'Update status to {status.value}'
        if analysis_type:
            job_name += f' (for {analysis_type})'

        j = b.new_job(job_name, job_attrs)
        j.image(image_path('cpg_workflows'))

        calc_size_script = None
        if outputs:
            if isinstance(outputs, dict):
                output_path = outputs[main_output_key]
            else:
                assert isinstance(outputs, str | Path)
                output_path = outputs
            assert isinstance(output_path, str | Path)
            calc_size_script = f"""
        from cloudpathlib import CloudPath
        meta['size'] = CloudPath('{str(output_path)}').stat().st_size
        """
        cmd = dedent(
            f"""\
        cat <<EOT >> update.py
        {calc_size_script}
        {self.updater_script(entry_id=entry_id, status=status)}
        EOT
        python3 update.py
        """
        )
        j.command(command(cmd, rm_leading_space=False, setup_gcp=True))
        return j

    @abstractmethod
    def updater_script(self, entry_id: int, status: AnalysisStatus) -> str:
        pass


class MetamistStateProvider(StateProvider):
    """
    Works through creating and updating Metamist's Analysis entries.
    """

    def __init__(self):
        super().__init__()

    @abstractmethod
    def read_state(
        self, cohort: Cohort, run_id: str
    ) -> dict[str, dict[str, AnalysisStatus]]:
        """
        On workflow-creating time, initialise state for each stage.
        Would read state for each stage+target into a dictionary, indexed by stage ID,
        then by target ID.
        """
        for dataset in cohort.get_datasets():
            gvcf_by_sid = get_metamist().get_analyses_by_sid(
                dataset.get_sample_ids(),
                analysis_type=AnalysisType.GVCF,
                dataset=dataset.name,
            )
            cram_by_sid = get_metamist().get_analyses_by_sid(
                dataset.get_sample_ids(),
                analysis_type=AnalysisType.CRAM,
                dataset=dataset.name,
            )
            for sample in dataset.get_samples():
                if (analysis := gvcf_by_sid.get(sample.id)) and analysis.output:
                    if analysis.output != sample.make_gvcf_path().path:
                        logging.warning(
                            f'GVCF path {analysis.output} does not match expected '
                            f'{sample.make_gvcf_path().path}'
                        )
                    sample.gvcf = analysis.output
                if (analysis := cram_by_sid.get(sample.id)) and analysis.output:
                    if analysis.output != sample.make_cram_path().path:
                        logging.warning(
                            f'CRAM path {analysis.output} does not match expected '
                            f'{sample.make_cram_path().path}'
                        )
                    sample.cram = analysis.output
        return {}

    def record_status(
        self,
        outputs: dict | str | Path | None,
        status: AnalysisStatus,
        stage_name: str,
        target: Target,
        analysis_type: str,
        meta: dict | None = None,
        dataset: str | None = None,
        main_output_key: str | None = None,
    ) -> int:
        """
        Record status as an Analysis entry
        """
        output_path: str | Path | None = None
        if isinstance(outputs, dict):
            output_path = outputs[main_output_key]
        elif isinstance(outputs, str | Path):
            output_path = outputs
        assert isinstance(output_path, str | Path | None)
        return get_metamist().create_analysis(
            output=outputs,
            type_=analysis_type,
            status=status,
            sample_ids=target.get_sample_ids(),
            meta=meta,
            dataset=dataset,
        )

    @abstractmethod
    def updater_script(self, entry_id: int, status: AnalysisStatus) -> str:
        return f"""
        from sample_metadata.apis import AnalysisApi
        from sample_metadata.models import AnalysisUpdateModel, AnalysisStatus
        from sample_metadata import exceptions
        import traceback

        meta = dict()
        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id={entry_id},
                analysis_update_model=AnalysisUpdateModel(
                    status=AnalysisStatus('{status.value}'),
                    meta=meta,
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        """


class FSStateProvider(StateProvider):
    """
    Works through checking stage's outputs existence on the file system.
    """


class JsonFileStateProvider(StateProvider):
    """
    Works through updating a JSON file.
    """
