"""
Metamist wrapper to report analysis progress.
"""
from abc import ABC, abstractmethod
from typing import Callable

from cpg_utils.config import get_config
from hailtop.batch.job import Job
from hailtop.batch import Batch

from .targets import Target


def complete_analysis_job(
        output: str,
        analysis_type: str,
        sg_ids: list[str],
        project_name: str,
        meta: dict,
        update_analysis_meta: Callable | None = None,
        tolerate_missing: bool = False
):
    """
    a job to be called within the batch as a pythonJob
    this will register the analysis outputs from a Stage

    Args:
        output (str): path to the output file
        analysis_type (str): metamist analysis type
        sg_ids (list[str]): all CPG IDs
        project_name (str): project/dataset name
        meta (dict): any metadata to add
        update_analysis_meta (Callable | None): function to update analysis meta
        tolerate_missing (bool): if True, allow missing output

    Returns:

    """
    import traceback
    from cpg_utils import to_path
    from metamist.apis import AnalysisApi
    from metamist.exceptions import ApiException
    from metamist.models import AnalysisStatus, Analysis
    assert isinstance(output, str)
    output_cloudpath = to_path(output)

    if update_analysis_meta is not None:
        meta | update_analysis_meta(output)

    # we know that es indexes are registered names, not files/dirs
    # skip all relevant checks for this output type
    if analysis_type != 'es-index':

        if not output_cloudpath.exists():
            if tolerate_missing:
                print(f"Output {output} doesn't exist, allowing silent return")
                return
            raise ValueError(f"Output {output} doesn't exist")

        # add file size to meta
        if not output_cloudpath.is_dir():
            meta |= {'size': output_cloudpath.stat().st_size}

    this_analysis = Analysis(
        type=analysis_type,
        status=AnalysisStatus('completed'),
        output=output,
        sequencing_group_ids=sg_ids,
        meta=meta,
    )
    aapi = AnalysisApi()
    try:
        a_id = aapi.create_analysis(project=project_name, analysis=this_analysis)
    except ApiException:
        traceback.print_exc()
        raise
    else:
        print(
            f'Created Analysis(id={a_id}, type={analysis_type}, '
            f'output={output}) in {project_name}'
        )
        return


class StatusReporterError(Exception):
    """
    Error thrown by StatusReporter.
    """


class StatusReporter(ABC):
    """
    Status reporter
    """

    @abstractmethod
    def create_analysis(
        self,
        b: Batch,
        output: str,
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        job_attr: dict | None = None,
        meta: dict | None = None,
        update_analysis_meta: Callable | None = None,
        tolerate_missing_output: bool = False,
        project_name: str | None = None,
    ):
        """
        Record analysis entry.
        """


class MetamistStatusReporter(StatusReporter):
    """
    Job status reporter. Works through creating metamist Analysis entries.
    """

    def __init__(self) -> None:
        super().__init__()

    def create_analysis(
        self,
        b: Batch,
        output: str,
        analysis_type: str,
        target: Target,
        jobs: list[Job] | None = None,
        job_attr: dict | None = None,
        meta: dict | None = None,
        update_analysis_meta: Callable | None = None,
        tolerate_missing_output: bool = False,
        project_name: str | None = None,
    ):
        """
        Create completed analysis job
        """

        # no jobs means no output, so no need to create analysis
        if not jobs:
            return

        if meta is None:
            meta = {}

        # find all relevant SG IDs
        sg_ids = target.get_sequencing_group_ids()
        py_job = b.new_python_job(
            f'Register analysis output {output}', job_attr or {} | {'tool': 'metamist'}
        )
        py_job.image(get_config()['workflow']['driver_image'])
        py_job.call(
            complete_analysis_job,
            str(output),
            analysis_type,
            sg_ids,
            project_name,
            meta,
            update_analysis_meta,
            tolerate_missing_output
        )

        py_job.depends_on(*jobs)
