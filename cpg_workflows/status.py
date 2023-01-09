"""
Workflow state provider.
"""
import inspect
from abc import ABC, abstractmethod
from typing import Callable, Any

from cpg_utils import Path
from cpg_utils.hail_batch import command, image_path
from hailtop.batch.job import Job
from hailtop.batch import Batch

from .targets import Target, Cohort
from .metamist import get_metamist, AnalysisStatus, AnalysisType


class StateProviderError(Exception):
    """
    Error thrown by StatusReporter.
    """


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
        update_analysis_meta: Callable[[str], dict] | None = None,
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

        meta_updaters_definitions = ''
        meta_updaters_funcs: list[Callable[[str], dict]] = []

        output_path = None
        if outputs:
            if isinstance(outputs, dict):
                output_path = outputs[main_output_key]
            else:
                assert isinstance(outputs, str | Path)
                output_path = outputs
            assert isinstance(output_path, str | Path)

            if isinstance(output_path, Path):
                meta_updaters_funcs.append(_calculate_size)
            if update_analysis_meta:
                meta_updaters_funcs.append(update_analysis_meta)

            for func in meta_updaters_funcs:
                definition = inspect.getsource(func)
                if not definition.startswith('def '):
                    raise StateProviderError(
                        f'Status updater must be a module-level function: {str(func)}'
                    )
                meta_updaters_definitions += definition + '\n'

        cmd = f"""
cat <<EOT >> update.py

from typing import Any, Callable

{meta_updaters_definitions}

{inspect.getsource(self.get_status_updater_function())}

{self.get_status_updater_function().__name__}(
    analysis_id={entry_id},
    new_status="{status.value}",
    updater_funcs=[{', '.join(f.__name__ for f in meta_updaters_funcs)}],
    output_path={'"' + str(output_path) + '"' if output_path else 'None'},
)

EOT
python3 update.py
"""

        j.command(command(cmd, rm_leading_space=False, setup_gcp=True))
        return j

    @abstractmethod
    def get_status_updater_function(self) -> Callable:
        pass


def _calculate_size(output_path: str) -> dict[str, Any]:
    """
    Self-contained function to calculate size of an object at given path.
    @param output_path: remote path of the output file
    @return: dictionary to merge into Analysis.meta
    """
    from cloudpathlib import CloudPath

    size = CloudPath(str(output_path)).stat().st_size
    return dict(size=size)


class MetamistStateProvider(StateProvider):
    """
    Job status reporter. Works through creating and updating metamist Analysis entries.
    """

    def __init__(self):
        super().__init__()

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
                    assert analysis.output == sample.make_gvcf_path().path, (
                        analysis.output,
                        sample.make_gvcf_path().path,
                    )
                    sample.gvcf = sample.make_gvcf_path()
                elif sample.make_gvcf_path().exists():
                    sample.gvcf = sample.make_gvcf_path()
                if (analysis := cram_by_sid.get(sample.id)) and analysis.output:
                    assert analysis.output == sample.make_cram_path().path, (
                        analysis.output,
                        sample.make_cram_path().path,
                    )
                    sample.cram = sample.make_cram_path()
                elif sample.make_cram_path().exists():
                    sample.cram = sample.make_cram_path()
        return {}

    def get_status_updater_function(self) -> Callable:
        """
        Get self-contained function that updates the status of a job.
        """

        def _update_analysis_status(
            analysis_id: int,
            new_status: str,
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
            from sample_metadata.model.analysis_status import (
                AnalysisStatus as MmAnalysisStatus,
            )
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
                        status=MmAnalysisStatus(new_status),
                        meta=meta,
                    ),
                )
            except exceptions.ApiException:
                traceback.print_exc()

        return _update_analysis_status

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


class JsonFileStateProvider(StateProvider):
    """
    Works through updating a JSON file.
    """
