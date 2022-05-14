"""
Helper functions to submit GATK-SV jobs.
"""

import json
import logging
import os
from os.path import join, dirname
from typing import Any

from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch import Batch

from cpg_pipes import images, Path
from cpg_pipes.hb.batch import make_job_name
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.targets import Dataset

from cpg_utils.hail_batch import reference_path, image_path

logger = logging.getLogger(__file__)


GATK_SV_COMMIT = 'd8467b2c2ac2eea9af07dcfb5cf3de46d719f7ef'
SV_CALLERS = ['manta', 'wham']

ACCESS_LEVEL = os.environ['CPG_ACCESS_LEVEL']


def get_dockers(keys: list[str]) -> dict[str, str]:
    """parse the WDL inputs with docker images"""
    with open(join(dirname(__file__), 'dockers.json')) as f:
        return {k: image_path(v) for k, v in json.load(f).items() if k in keys}


def get_references(keys: list[str]) -> dict[str, str]:
    """parse the WDL inputs with reference files"""
    with open(join(dirname(__file__), 'references.json')) as f:
        return {k: reference_path(v) for k, v in json.load(f).items() if k in keys}


def add_gatksv_job(
    batch: Batch,
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any], 
    expected_out_dict: dict[str, Path],
    sample_id: str | None = None,
):
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output. 
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sample_id: 
        output_prefix = join(output_prefix, sample_id)

    outputs_to_collect = dict()
    for key in expected_out_dict.keys():
        outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset.name)
    assert ACCESS_LEVEL
    output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=dataset.stack,
        access_level=ACCESS_LEVEL,
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict={f'{wfl_name}.{k}': v for k, v in input_dict.items()},
        outputs_to_collect=outputs_to_collect,
        driver_image=images.DRIVER_IMAGE,
    )

    copy_j = batch.new_job(f'{job_prefix}: copy outputs')
    copy_j.image(images.DRIVER_IMAGE)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        cmds.append(f'gsutil cp $(cat {resource}) {out_path}')
    copy_j.command(wrap_command(cmds, setup_gcp=True))

    return output_dict, copy_j
