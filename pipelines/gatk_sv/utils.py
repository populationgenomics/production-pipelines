"""
Helper functions to submit GATK-SV jobs.
"""

import json
import logging
import os
from os.path import join, dirname
from typing import Any, TypedDict

from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch import Batch

from cpg_pipes import images, Path, to_path
from cpg_pipes.hb.batch import make_job_name
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.targets import Dataset

from cpg_utils.hail_batch import reference_path, image_path

logger = logging.getLogger(__file__)


GATK_SV_COMMIT = 'e4c5ffb1596e0de93eb491ae6d7cb6ec46874d11'
SV_CALLERS = ['manta', 'wham']


def get_dockers(keys: list[str]) -> dict[str, str]:
    """
    Dict of WDL inputs with docker image paths
    """
    with open(join(dirname(__file__), 'dockers.json')) as f:
        return {k: image_path(v) for k, v in json.load(f).items() if k in keys}


def _to_path(val: str) -> str:
    """
    Relative paths are resolved with respect to CPG_REFERENCE_PATH
    """
    if to_path(val).is_absolute():
        return val
    return str(reference_path(val))


def get_gcnv_models() -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with gCNV models
    """
    res: dict[str, str | list[str]] = {}
    with open(join(dirname(__file__), 'gcnv.json')) as f:
        d = json.load(f)
        res['contig_ploidy_model_tar'] = _to_path(d['contig_ploidy_model_tar'])
        res['gcnv_model_tars'] = [
            _to_path(d['model_tar_tmpl'].format(shard=i)) 
            for i in range(d['model_tar_cnt'])
        ]
        res['ref_panel_samples'] = d['ref_panel_samples']
        res['ref_panel_PE_files'] = [
            _to_path(d['ref_panel_PE_file_tmpl'].format(sample=s))
            for s in d['ref_panel_samples']
        ]
        res['ref_panel_SE_files'] = [
            _to_path(d['ref_panel_SE_file_tmpl'].format(sample=s))
            for s in d['ref_panel_samples']
        ]
    return res


def get_references(
    keys: list[str | dict[str, str]] | None = None
) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """    
    def _val_to_ref_path(val: str | list[str]) -> str | list[str]:
        """Expecting a path or a list of paths"""
        if isinstance(val, list):
            return [_to_path(x) for x in val]
        else:
            return _to_path(val)

    with open(join(dirname(__file__), 'references.json')) as f:
        ref_d = json.load(f)
    inputs_d = {}
    if keys:
        for key in keys:
            # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
            if isinstance(key, dict):
                key, ref_d_key = list(key.items())[0]
            else:
                ref_d_key = key
            # e.g. GATKSVPipelineBatch.rmsk -> rmsk
            ref_d_key = ref_d_key.split('.')[-1] 
            inputs_d[key] = ref_d[ref_d_key]
    else:
        inputs_d = ref_d
    inputs_d = {k: _val_to_ref_path(v) for k, v in inputs_d.items()}
    return inputs_d


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
    access_level = os.environ['CPG_ACCESS_LEVEL']
    assert access_level
    output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=dataset.stack,
        access_level=access_level,
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
