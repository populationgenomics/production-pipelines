"""
Workflow state provider.
"""
import inspect
import json
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
        update_analysis_meta: Callable[[str], dict] | None = None,
    ) -> list[Job]:
        """
        Record QUEUED status for a stage, and insert jobs that update status to
        IN_PROGRESS and COMPLETED.
        """
        if not jobs:
            return []

        self.record_status(
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
            stage_name=stage_name,
            target=target,
            status=AnalysisStatus.IN_PROGRESS,
            job_attrs=(job_attrs or {}) | {'tool': 'update_state'},
        )
        # 2. Queue a job that updates the status to COMPLETED
        completed_j = self.add_status_updater_job(
            b,
            stage_name=stage_name,
            target=target,
            status=AnalysisStatus.COMPLETED,
            job_attrs=(job_attrs or {}) | {'tool': 'update_state'},
            outputs=outputs,
            main_output_key=main_output_key,
            update_analysis_meta=update_analysis_meta,
        )

        in_progress_j.depends_on(*(prev_jobs or []))
        completed_j.depends_on(*jobs)
        return [in_progress_j, *jobs, completed_j]

    @abstractmethod
    def get_entry_id(self, stage_name: str, target: Target) -> str:
        pass

    def add_status_updater_job(
        self,
        b: Batch,
        stage_name: str,
        target: Target,
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

        entry_id = self.get_entry_id(stage_name, target)

        cmd = f"""
cat <<EOT >> update.py

from typing import Any, Callable
from cpg_workflows.status import {self.__class__.__name__}

{meta_updaters_definitions}

{self.__class__.__name__}.update_status(
    entry_id="{entry_id}",
    new_status="{status.value}",
    updater_funcs=[{', '.join(f.__name__ for f in meta_updaters_funcs)}],
    output_path={'"' + str(output_path) + '"' if output_path else 'None'},
)

EOT
python3 update.py
"""

        j.command(command(cmd, rm_leading_space=False, setup_gcp=True))
        return j

    @staticmethod
    @abstractmethod
    def update_status(
        entry_id: str,
        new_status: str,
        updater_funcs: list[Callable[[str], dict[str, Any]]] | None = None,
        output_path: str | None = None,
    ):
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
        self.analysis_entry_ids: dict[tuple[str, str], str] = {}

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

    def get_entry_id(self, stage_name: str, target: Target) -> str:
        """
        Get entry ID for a given stage and target.
        """
        entry_id = self.analysis_entry_ids.get((stage_name, target.target_id))
        if not entry_id:
            raise StateProviderError(f'No entry ID for {stage_name} {target.target_id}')
        try:
            int_entry_id = int(entry_id)
        except ValueError:
            raise StateProviderError(
                f'Metamist entry ID must be numerical, got {entry_id}'
            )
        return entry_id

    @staticmethod
    def update_status(
        entry_id: str,
        new_status: str,
        updater_funcs: list[Callable[[str], dict[str, Any]]] | None = None,
        output_path: str | None = None,
    ) -> None:
        """
        Self-contained function to update Metamist analysis entry.
        @param entry_id: ID of the Analysis entry
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
                analysis_id=int(entry_id),
                analysis_update_model=AnalysisUpdateModel(
                    status=MmAnalysisStatus(new_status),
                    meta=meta,
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()

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
    ):
        """
        Record status as an Analysis entry
        """
        output_path: str | Path | None = None
        if isinstance(outputs, dict):
            if not main_output_key:
                raise StateProviderError(
                    f'Cannot create Analysis: `analysis_key` '
                    f'must be set with the @stage decorator to select value from '
                    f'the expected_outputs dict: {outputs}'
                )
            if main_output_key not in outputs:
                raise StateProviderError(
                    f'Cannot create Analysis: `analysis_key` '
                    f'"{main_output_key}" is not found in the expected_outputs '
                    f'dict {outputs}'
                )
            output_path = outputs[main_output_key]

        assert isinstance(output_path, str | Path | None)
        analysis_id = get_metamist().create_analysis(
            output=outputs,
            type_=analysis_type,
            status=status,
            sample_ids=target.get_sample_ids(),
            meta=meta,
            dataset=dataset,
        )
        self.analysis_entry_ids[(stage_name, target.target_id)] = str(analysis_id)


class JsonFileStateProvider(StateProvider):
    """
    Works through updating a JSON file.
    """

    def __init__(self, prefix: Path):
        super().__init__()
        self.prefix = prefix

    def read_state(
        self, cohort: Cohort, run_id: str
    ) -> dict[str, dict[str, AnalysisStatus]]:
        pass

    def make_json_path(self, stage_name: str, target: Target) -> Path:
        return self.prefix / f'{stage_name}-{target.target_id}.json'

    def get_entry_id(self, stage_name: str, target: Target) -> str:
        return str(self.make_json_path(stage_name, target))

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
        if isinstance(outputs, dict):
            outputs = {k: str(v) for k, v in outputs.items()}
        else:
            outputs = str(outputs)
        d = dict(
            outputs=outputs,
            type=analysis_type,
            status=status.value,
            sample_ids=target.get_sample_ids(),
            meta=meta or {},
        )
        if dataset:
            d['dataset'] = dataset

        path = self.make_json_path(stage_name, target)
        with path.open('w') as f:
            json.dump(d, f)
            return 0

    @staticmethod
    def update_status(
        entry_id: str,
        new_status: str,
        updater_funcs: list[Callable[[str], dict[str, Any]]] | None = None,
        output_path: str | None = None,
    ) -> None:
        """
        Self-contained function to update Metamist analysis entry.
        @param entry_id: ID of the status entry
        @param new_status: new status to assign to the entry
        @param updater_funcs: list of functions to update the entry's metadata,
        assuming output_path as input parameter
        @param output_path: remote path of the output file, to be passed to the updaters
        """
        import json
        from cpg_utils import to_path

        path = to_path(entry_id)

        meta: dict[str, Any] = dict()
        if output_path and updater_funcs:
            for func in updater_funcs or []:
                meta |= func(output_path)

        with path.open('r') as f:
            d = json.load(f)
        d['status'] = new_status
        d['meta'] |= meta
        with path.open('w') as f:
            json.dump(d, f)
