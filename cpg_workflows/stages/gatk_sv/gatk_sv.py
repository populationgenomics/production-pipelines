"""
Driver for computing structural variants from GATK-SV from WGS data.
"""
from os.path import join
from typing import Any

from analysis_runner.cromwell import (
    run_cromwell_workflow_from_repo_and_get_outputs,
    CromwellOutputType,
)
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, ConfigError
from cpg_utils.hail_batch import command, reference_path, image_path
from cpg_workflows.batch import make_job_name, Batch, get_batch
from cpg_workflows.workflow import (
    stage,
    StageOutput,
    DatasetStage,
    StageInput,
    Dataset,
    Cohort,
    CohortStage
)

GATK_SV_COMMIT = 'a73237cf9d9e321df3aa81c890def7b504a25c7f'
SV_CALLERS = ['manta', 'wham', 'scramble']
_FASTA = None


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
    sample_id: str | None = None,
    driver_image: str | None = None,
) -> list[Job]:
    """
    Generic function to add a job that would run one GATK-SV workflow.
    """
    # Where Cromwell writes the output.
    # Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/{wfl_name}/{dataset.name}'
    if sample_id:
        output_prefix = join(output_prefix, sample_id)

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

    job_prefix = make_job_name(wfl_name, sample=sample_id, dataset=dataset.name)
    submit_j, output_dict = run_cromwell_workflow_from_repo_and_get_outputs(
        b=batch,
        job_prefix=job_prefix,
        dataset=get_config()['workflow']['dataset'],
        access_level=get_config()['workflow']['access_level'],
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow=f'{wfl_name}.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=driver_image,
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
                reference_path('gatk_sv/ref_panel_bincov_matrix')
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
    with combined_ped_path.open('w') as out:
        with cohort.write_ped_file().open() as f:
            out.write(f.read())
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with reference_path('gatk_sv/ped_file').open() as f:
            out.write(f.read())
    return combined_ped_path


# @stage(required_stages=GenotypeBatch)
# class RegenotypeCNVs(DatasetStage):
#     """
#     Optional, Re-genotypes probable mosaic variants across multiple batches.
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict:
#         pass
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
#         pass


@stage(required_stages=[FilterBatch, GenotypeBatch, GatherBatchEvidence])
class MakeCohortVcf(DatasetStage):
    """
    Combines variants across multiple batches, resolves complex variants, re-genotypes,
    and performs final VCF clean-up.

    If RegenotypeCNVs is run, this stage will use the output of that stage as input.
    Otherwise, it will use the output of GenotypeBatch.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        """

        Args:
            dataset ():

        Returns:

        """
        ending_by_key = {
            'vcf': '.cleaned.vcf.gz',
            'vcf_index': '.cleaned.vcf.gz.tbi',
            'vcf_qc': '.cleaned_SV_VCF_QC_output.tar.gz',
            # if merge_intermediate_vcfs is enabled
            # 'cluster_vcf': '.combine_batches.vcf.gz',
            # 'cluster_vcf_index': '.combine_batches.vcf.gz.tbi',
            # 'complex_resolve_vcf': '.complex_resolve.vcf.gz',
            # 'complex_resolve_vcf_index': '.complex_resolve.vcf.gz.tbi',
            # 'complex_genotype_vcf': '.complex_genotype.vcf.gz',
            # 'complex_genotype_vcf_index': '.complex_genotype.vcf.gz.tbi',
            'metrics_file_makecohortvcf': '.metrics.tsv'
        }
        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """

        Args:
            dataset ():
            inputs ():

        Returns:

        """
        batchevidence_d = inputs.as_dict(dataset, GatherBatchEvidence)
        filterbatch_d = inputs.as_dict(dataset, FilterBatch)
        genotypebatch_d = inputs.as_dict(dataset, GenotypeBatch)

        input_dict: dict[str, Any] = {
            'cohort_name': dataset.name,
            'batches': [dataset.name],
            'ped_file': make_combined_ped(dataset),
            'ref_dict': str(get_fasta().with_suffix('.dict')),
            'chr_x': 'chrX',
            'chr_y': 'chrY',
            'min_sr_background_fail_batches': 0.5,
            'max_shard_size_resolve': 500,
            'max_shards_per_chrom_clean_vcf_step1': 200,
            'min_records_per_shard_clean_vcf_step1': 5000,
            'clean_vcf1b_records_per_shard': 10000,
            'samples_per_clean_vcf_step2_shard': 100,
            'clean_vcf5_records_per_shard': 5000,
            'random_seed': 0,
            # not explicit, but these VCFs require indices
            'pesr_vcfs': [genotypebatch_d['genotyped_pesr_vcf']],
            'depth_vcfs': [genotypebatch_d['genotyped_depth_vcf']],
            'disc_files': [batchevidence_d['merged_PE']],
            'bincov_files': [batchevidence_d['merged_bincov']],
            'raw_sr_bothside_pass_files': [genotypebatch_d['sr_bothside_pass']],
            'raw_sr_background_fail_files': [genotypebatch_d['sr_background_fail']],
            'depth_gt_rd_sep_files': [genotypebatch_d['trained_genotype_depth_depth_sepcutoff']],
            'median_coverage_files': [batchevidence_d['median_cov']],
            'rf_cutoff_files': [filterbatch_d['cutoffs']],
        }

        input_dict |= get_references(
            [
                'bin_exclude',
                'mei_bed',
                'depth_exclude_list',
                'empty_file',
                # same attr, two names
                'primary_contigs_list',
                {'contig_list': 'primary_contigs_list'},
                {'allosome_fai': 'allosome_file'},
                {'cytobands': 'cytoband'},
                {'pe_exclude_list': 'pesr_exclude_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'linux_docker',
            ]
        )
        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)


@stage(required_stages=MakeCohortVcf)
class AnnotateVcf(DatasetStage):
    """
    Add annotations, such as the inferred function and allele frequencies of variants,
    to final VCF.

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
    """

    def expected_outputs(self, dataset: Dataset) -> dict:
        # vcf_path = str(self.tmp_prefix / 'sv' / f'{dataset.name}.vcf.gz')
        ending_by_key = {
            'output_vcf': '.annotated.vcf.gz',
            'output_vcf_idx': '.annotated.vcf.gz.tbi',
        }
        d: dict[str, Path] = {}
        for key, ending in ending_by_key.items():
            fname = f'{dataset.name}{ending}'
            d[key] = dataset.prefix() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:

        make_vcf_d = inputs.as_dict(dataset, MakeCohortVcf)

        input_dict: dict[str, Any] = {
            'prefix': dataset.name,
            'vcf': make_vcf_d['vcf'],
            'vcf_idx': make_vcf_d['vcf_index'],
            'ped_file': make_combined_ped(dataset),
            'sv_per_shard': 5000,
            'max_shards_per_chrom_step1': 200,
            'min_records_per_shard_step1': 5000
        }
        input_dict |= get_references(
            [
                'protein_coding_gtf',
                {'contig_list': 'primary_contigs_list'},
            ]
        )

        # images!
        input_dict |= get_images(
            [
                'sv_pipeline_docker',
                'sv_base_mini_docker',
                'gatk_docker',
            ]
        )
        expected_d = self.expected_outputs(dataset)
        jobs = add_gatk_sv_jobs(
            batch=get_batch(),
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )
        return self.make_outputs(dataset, data=expected_d, jobs=jobs)
