"""
Common methods for all GATK-SV workflows
"""

import re
from enum import Enum
from functools import lru_cache
from os.path import join
from random import randint
from typing import Any

from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import ConfigError, config_retrieve, image_path, reference_path
from cpg_utils.cromwell import CromwellOutputType, run_cromwell_workflow_from_repo_and_get_outputs
from cpg_utils.hail_batch import command, get_batch
from cpg_workflows.batch import make_job_name
from cpg_workflows.workflow import Cohort, Dataset, MultiCohort

GATK_SV_COMMIT = 'dc145a52f76a6f425ac3f481171040e78c0cfeea'
SV_CALLERS = ['manta', 'wham', 'scramble']
_FASTA = None
PED_FAMILY_ID_REGEX = re.compile(r'(^[A-Za-z0-9_]+$)')


class CromwellJobSizes(Enum):
    """
    Enum for polling intervals
    """

    SMALL = 'small'
    MEDIUM = 'medium'
    LARGE = 'large'


@lru_cache(maxsize=1)
def create_polling_intervals() -> dict:
    """
    Set polling intervals for cromwell status
    these values are integers, indicating seconds
    for each job size, there is a min and max value
    analysis-runner implements a backoff-retrier when checking for
    success, with a minimum value, gradually reaching a max ceiling

    a config section containing overrides would look like

    [cromwell_polling_intervals.medium]
    min = 69
    max = 420
    """

    # create this dict with default values
    polling_interval_dict = {
        CromwellJobSizes.SMALL: {'min': 30, 'max': 140},
        CromwellJobSizes.MEDIUM: {'min': 40, 'max': 400},
        CromwellJobSizes.LARGE: {'min': 200, 'max': 2000},
    }

    # update if these exist in config
    for job_size in CromwellJobSizes:
        if val := config_retrieve(['cromwell_polling_intervals', job_size.value], False):
            polling_interval_dict[job_size].update(val)
    return polling_interval_dict


def get_fasta() -> Path:
    """
    find or return the fasta to use
    """
    global _FASTA
    if _FASTA is None:
        fasta_path = config_retrieve(['workflow', 'ref_fasta'], False) or reference_path('broad/ref_fasta')
        _FASTA = to_path(fasta_path)
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
    image_keys = config_retrieve(['images']).keys()

    if not allow_missing:
        query_keys = set(keys)
        if not query_keys.issubset(image_keys):
            raise KeyError(f'Unknown image keys: {query_keys - image_keys}')

    return {k: image_path(k) for k in image_keys if k in keys}


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
            res[key] = reference_path(f'gatk_sv/{ref_d_key}')
        except KeyError:
            res[key] = reference_path(f'broad/{ref_d_key}')
        except ConfigError:
            res[key] = reference_path(f'broad/{ref_d_key}')

    return res


def add_gatk_sv_jobs(
    dataset: Dataset,
    wfl_name: str,
    # "dict" is invariant (supports updating), "Mapping" is covariant (read-only)
    # we have to support inputs of type dict[str, str], so using Mapping here:
    input_dict: dict[str, Any],
    expected_out_dict: dict[str, Path | list[Path]],
    sequencing_group_id: str | None = None,
    driver_image: str | None = None,
    labels: dict[str, str] | None = None,
    job_size: CromwellJobSizes = CromwellJobSizes.MEDIUM,
) -> list[Job]:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """

    # create/retrieve dictionary of polling intervals for cromwell status
    polling_intervals = create_polling_intervals()

    # obtain upper and lower polling bounds for this job size
    polling_minimum = randint(polling_intervals[job_size]['min'], polling_intervals[job_size]['min'] * 2)
    polling_maximum = randint(polling_intervals[job_size]['max'], polling_intervals[job_size]['max'] * 2)

    # If a config section exists for this workflow, apply overrides
    if override := config_retrieve(['resource_overrides', wfl_name], False):
        input_dict |= override

    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sequencing_group_id:
        output_prefix = join(output_prefix, sequencing_group_id)

    outputs_to_collect = dict()
    for key, value in expected_out_dict.items():
        if isinstance(value, list):
            outputs_to_collect[key] = CromwellOutputType.array_path(name=f'{wfl_name}.{key}', length=len(value))
        else:
            outputs_to_collect[key] = CromwellOutputType.single_path(f'{wfl_name}.{key}')

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

    job_prefix = make_job_name(wfl_name, sequencing_group=sequencing_group_id, dataset=dataset.name)

    submit_j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=get_batch(),
        job_prefix=job_prefix,
        dataset=config_retrieve(['workflow', 'dataset']),
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
        copy_outputs_to_gcp=config_retrieve(['workflow', 'copy_outputs'], False),
        labels=labels,
        min_watch_poll_interval=polling_minimum,
        max_watch_poll_interval=polling_maximum,
        time_limit_seconds=config_retrieve(['workflow', 'time_limit_seconds'], None),
    )

    copy_j = get_batch().new_job(f'{job_prefix}: copy outputs')
    copy_j.image(driver_image)
    cmds = []
    for key, resource in output_dict.items():
        out_path = expected_out_dict[key]
        if isinstance(resource, list):
            for source, dest in zip(resource, out_path):
                cmds.append(f'gcloud storage cp "$(cat {source})" "{dest}"')
        else:
            cmds.append(f'gcloud storage cp "$(cat {resource})" "{out_path}"')
    copy_j.command(command(cmds, setup_gcp=True))
    return [submit_j, copy_j]


def get_ref_panel(keys: list[str] | None = None) -> dict:
    # mandatory config entry
    ref_panel_samples = config_retrieve(['sv_ref_panel', 'ref_panel_samples'])
    return {
        k: v
        for k, v in {
            'ref_panel_samples': ref_panel_samples,
            'ref_panel_bincov_matrix': reference_path('broad/ref_panel_bincov_matrix'),
            'contig_ploidy_model_tar': reference_path('gatk_sv/contig_ploidy_model_tar'),
            'gcnv_model_tars': [
                str(reference_path('gatk_sv/model_tar_tmpl')).format(shard=i)
                for i in range(config_retrieve(['sv_ref_panel', 'model_tar_cnt']))
            ],
            'ref_panel_PE_files': [
                reference_path('gatk_sv/ref_panel_PE_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
            'ref_panel_SR_files': [
                reference_path('gatk_sv/ref_panel_SR_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
            'ref_panel_SD_files': [
                reference_path('gatk_sv/ref_panel_SD_file_tmpl').format(sample=s) for s in ref_panel_samples
            ],
        }.items()
        if not keys or k in keys
    }


def clean_ped_family_ids(ped_line: str) -> str:
    """
    Takes each line in the pedigree and cleans it up
    If the family ID already conforms to expectations, no action
    If the family ID fails, replace all non-alphanumeric/non-underscore
    characters with underscores

    >>> clean_ped_family_ids('family1\tchild1\t0\t0\t1\t0\\n')
    'family1\tchild1\t0\t0\t1\t0\\n'
    >>> clean_ped_family_ids('family-1-dirty\tchild1\t0\t0\t1\t0\\n')
    'family_1_dirty\tchild1\t0\t0\t1\t0\\n'

    Args:
        ped_line (str): line from the pedigree file, unsplit

    Returns:
        the same line with a transformed family id
    """

    split_line = ped_line.rstrip().split('\t')

    if re.match(PED_FAMILY_ID_REGEX, split_line[0]):
        return ped_line

    # if the family id is not valid, replace failing characters with underscores
    split_line[0] = re.sub(r'[^A-Za-z0-9_]', '_', split_line[0])

    # return the rebuilt string, with a newline at the end
    return '\t'.join(split_line) + '\n'


def make_combined_ped(cohort: Cohort | MultiCohort, prefix: Path) -> Path:
    """
    Create cohort + ref panel PED.
    Concatenating all samples across all datasets with ref panel

    See #578 - there are restrictions on valid characters in PED file
    """
    combined_ped_path = prefix / 'ped_with_ref_panel.ped'
    conf_ped_path = get_references(['ped_file'])['ped_file']
    assert isinstance(conf_ped_path, str)
    with combined_ped_path.open('w') as out:
        with cohort.write_ped_file().open() as f:
            # layer of family ID cleaning
            for line in f:
                out.write(clean_ped_family_ids(line))
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with to_path(conf_ped_path).open() as f:
            out.write(f.read())
    return combined_ped_path


def queue_annotate_sv_jobs(
    multicohort: MultiCohort,
    prefix: Path,
    input_vcf: Path,
    outputs: dict,
    labels: dict[str, str] | None = None,
) -> list[Job] | Job | None:
    """
    Helper function to queue jobs for SV annotation
    Enables common access to the same Annotation WDL for CNV & SV
    """
    input_dict: dict[str, Any] = {
        'vcf': input_vcf,
        'prefix': multicohort.name,
        'ped_file': make_combined_ped(multicohort, prefix),
        'sv_per_shard': 5000,
        'external_af_population': config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'external_af_ref_prefix': config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config_retrieve(['references', 'gnomad_sv']),
        'use_hail': False,
    }

    input_dict |= get_references(
        [
            'noncoding_bed',
            'protein_coding_gtf',
            {'contig_list': 'primary_contigs_list'},
        ],
    )

    # images!
    input_dict |= get_images(['sv_pipeline_docker', 'sv_base_mini_docker', 'gatk_docker'])
    jobs = add_gatk_sv_jobs(
        dataset=multicohort.analysis_dataset,
        wfl_name='AnnotateVcf',
        input_dict=input_dict,
        expected_out_dict=outputs,
        labels=labels,
    )
    return jobs


def queue_annotate_strvctvre_job(
    input_vcf,
    output_path: str,
    job_attrs: dict | None = None,
    name: str = 'AnnotateVcfWithStrvctvre',
) -> Job:
    """

    Args:
        input_vcf (ResourceFile): part of a resource group with the corresponding index
        output_path ():
        job_attrs (dict|None): job attributes
        name (str): name of the job

    Returns:
        The Strvctvre job
    """

    job_attrs = job_attrs or {}
    strv_job = get_batch().new_job('StrVCTVRE', job_attrs | {'tool': 'strvctvre'})

    strv_job.image(image_path('strvctvre'))
    strv_job.storage(config_retrieve(['resource_overrides', name, 'storage'], '10Gi'))
    strv_job.memory(config_retrieve(['resource_overrides', name, 'memory'], '16Gi'))

    strvctvre_phylop = get_references(['strvctvre_phylop'])['strvctvre_phylop']
    assert isinstance(strvctvre_phylop, str)

    local_phylop = get_batch().read_input(strvctvre_phylop)

    strv_job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})

    # run strvctvre
    strv_job.command(
        f'python StrVCTVRE.py -i {input_vcf} -o {strv_job.output["vcf.gz"]} -f vcf -p {local_phylop}',  # type: ignore
    )
    strv_job.command(f'tabix {strv_job.output["vcf.gz"]}')  # type: ignore

    get_batch().write_output(strv_job.output, str(output_path).replace('.vcf.gz', ''))
    return strv_job
