"""
Common methods for all GATK-SV workflows
"""

from os.path import join
from typing import Any

from analysis_runner.cromwell import (
    CromwellOutputType,
    run_cromwell_workflow_from_repo_and_get_outputs,
)
from cpg_utils import Path, to_path
from cpg_utils.config import ConfigError, get_config
from cpg_utils.hail_batch import Batch, command, image_path, reference_path
from hailtop.batch.job import Job

from cpg_workflows.batch import make_job_name
from cpg_workflows.workflow import Cohort, Dataset

GATK_SV_COMMIT = '6d6100082297898222dfb69fcf941d373d78eede'
SV_CALLERS = ['manta', 'wham', 'scramble']
_FASTA = None


def _sv_batch_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Callable, add meta[type] to custom analysis object
    """
    return {'type': 'gatk-sv-batch-calls'}


def _sv_filtered_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Callable, add meta[type] to custom analysis object
    """
    return {'type': 'gatk-sv-filtered-calls'}


def _sv_individual_meta(
    output_path: str,  # pylint: disable=W0613:unused-argument
) -> dict[str, Any]:
    """
    Callable, add meta[type] to custom analysis object
    """
    return {'type': 'gatk-sv-sequence-group-calls'}


def get_fasta() -> Path:
    """
    find or return the fasta to use
    """
    global _FASTA
    if _FASTA is None:
        _FASTA = to_path(
            get_config()['workflow'].get('ref_fasta')
            or reference_path('broad/ref_fasta')
        )
    return _FASTA


def get_images(keys: list[str], allow_missing=False) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths.

    Args:
        keys (list): all the images to get
        allow_missing (bool): if False, require all query keys to be found

    Returns:
        dict of image keys to image paths
        or AssertionError
    """

    if not allow_missing:
        image_keys = get_config()['images'].keys()
        query_keys = set(keys)
        if not query_keys.issubset(image_keys):
            raise KeyError(f'Unknown image keys: {query_keys - image_keys}')

    return {k: image_path(k) for k in get_config()['images'].keys() if k in keys}


def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        if isinstance(key, dict):
            key, ref_d_key = list(key.items())[0]
        else:
            ref_d_key = key
        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = str(reference_path(f'gatk_sv/{ref_d_key}'))
        except KeyError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))
        except ConfigError:
            res[key] = str(reference_path(f'broad/{ref_d_key}'))

    return res


def add_gatk_sv_jobs(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path | list[Path]],
    sequencing_group_id: str | None = None,
    driver_image: str | None = None,
    labels: dict[str, str] | None = None,
    cromwell_status_poll_interval: int = 60,
) -> list[Job]:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sequencing_group_id:
        output_prefix = join(output_prefix, sequencing_group_id)

    outputs_to_collect = dict()
    for key, value in expected_out_dict.items():
        if isinstance(value, list):
            outputs_to_collect[key] = CromwellOutputType.array_path(
                name=f'{wfl_name}.{key}', length=len(value)
            )
        else:
            outputs_to_collect[key] = CromwellOutputType.single_path(
                f'{wfl_name}.{key}'
            )

    driver_image = driver_image or image_path('cpg_workflows')

    # pre-process input_dict
    paths_as_strings: dict = {}
    for key, value in input_dict.items():
        if isinstance(value, Path):
            paths_as_strings[f'{wfl_name}.{key}'] = str(value)
        elif isinstance(value, (list, set)):
            paths_as_strings[f'{wfl_name}.{key}'] = [str(v) for v in value]
        else:
            paths_as_strings[f'{wfl_name}.{key}'] = value

    job_prefix = make_job_name(
        wfl_name, sequencing_group=sequencing_group_id, dataset=dataset.name
    )

    # config toggle decides if outputs are copied out
    copy_outputs = get_config()['workflow'].get('copy_outputs', False)

    submit_j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=get_config()['workflow']['dataset'],
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
        copy_outputs_to_gcp=copy_outputs,
        labels=labels,
        max_watch_poll_interval=cromwell_status_poll_interval,
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(driver_image)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        if isinstance(resource, list):
            for source, dest in zip(resource, out_path):
                cmds.append(f'gsutil cp "$(cat {source})" "{dest}"')
        else:
            cmds.append(f'gsutil cp "$(cat {resource})" "{out_path}"')
    copy_j.command(command(cmds, setup_gcp=True))
    return [submit_j, copy_j]


def get_ref_panel(keys: list[str] | None = None) -> dict:
    return {
        k: v
        for k, v in {
            'ref_panel_samples': get_config()['sv_ref_panel']['ref_panel_samples'],
            'ref_panel_bincov_matrix': str(
                reference_path('broad/ref_panel_bincov_matrix')
            ),
            'contig_ploidy_model_tar': str(
                reference_path('gatk_sv/contig_ploidy_model_tar')
            ),
            'gcnv_model_tars': [
                str(reference_path('gatk_sv/model_tar_tmpl')).format(shard=i)
                for i in range(get_config()['sv_ref_panel']['model_tar_cnt'])
            ],
            'ref_panel_PE_files': [
                str(reference_path('gatk_sv/ref_panel_PE_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SR_files': [
                str(reference_path('gatk_sv/ref_panel_SR_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
            'ref_panel_SD_files': [
                str(reference_path('gatk_sv/ref_panel_SD_file_tmpl')).format(sample=s)
                for s in get_config()['sv_ref_panel']['ref_panel_samples']
            ],
        }.items()
        if not keys or k in keys
    }


def make_combined_ped(cohort: Cohort, prefix: Path) -> Path:
    """
    Create cohort + ref panel PED.
    Concatenating all samples across all datasets with ref panel
    """
    combined_ped_path = prefix / 'ped_with_ref_panel.ped'
    conf_ped_path = get_references(['ped_file'])['ped_file']
    with combined_ped_path.open('w') as out:
        with cohort.write_ped_file().open() as f:
            out.write(f.read())
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with to_path(conf_ped_path).open() as f:
            out.write(f.read())
    return combined_ped_path
