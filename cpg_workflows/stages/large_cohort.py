import json
from logging import config
from typing import TYPE_CHECKING, Any, Final, Tuple

from numpy import require
from requests import get

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, get_config, image_path, output_path
from cpg_utils.hail_batch import get_batch, query_command
from cpg_workflows.jobs.vep import add_vep_jobs
from cpg_workflows.large_cohort.combiner import combiner
from cpg_workflows.targets import Cohort
from cpg_workflows.utils import get_all_fragments_from_manifest, slugify
from cpg_workflows.workflow import (
    CohortStage,
    StageInput,
    StageOutput,
    get_workflow,
    stage,
)

if TYPE_CHECKING:
    from graphql import DocumentNode

    from hailtop.batch.job import PythonJob

HAIL_QUERY: Final = 'hail query'

SHARD_MANIFEST = 'shard-manifest.txt'


# TODO, update analysis_meta here to pull the gvcf.type, and store this in metamist.
@stage(analysis_type='combiner', analysis_keys=['vds'])
class Combiner(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Any]:
        combiner_config: dict[str, str] = config_retrieve('combiner')

        # Allow user to specify a custom VDS path
        vds_path = combiner_config.get('vds_path', False)
        if not vds_path:
            output_vds_name: str = slugify(
                f'{cohort.id}-{combiner_config["vds_version"]}',
            )
            vds_path = cohort.analysis_dataset.prefix() / 'vds' / f'{cohort.name}' / f'{output_vds_name}.vds'
        else:
            vds_path = to_path(vds_path)

        return {'vds': vds_path}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        output_paths = self.expected_outputs(cohort)

        # include the list of all VDS IDs in the plan name
        if vds_ids := config_retrieve(['combiner', 'vds_analysis_ids']):
            ids_list_as_string: str = '_'.join([str(id) for id in sorted(vds_ids)])
            combiner_plan_name: str = f'combiner_{ids_list_as_string}'
        else:
            combiner_plan_name = f'combiner-{cohort.name}'

        combiner_plan: str = str(self.get_stage_cohort_prefix(cohort, 'tmp') / f'{combiner_plan_name}.json')

        j: PythonJob = combiner(
            cohort=cohort,
            output_vds_path=str(output_paths['vds']),
            save_path=combiner_plan,
        )

        return self.make_outputs(cohort, output_paths, [j])


@stage(required_stages=[Combiner])
class SampleQC(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if sample_qc_version := config_retrieve(['large_cohort', 'output_versions', 'sample_qc'], default=None):
            sample_qc_version = slugify(sample_qc_version)

        sample_qc_version = sample_qc_version or get_workflow().output_version
        sample_qc_path = cohort.analysis_dataset.prefix() / get_workflow().name / sample_qc_version / 'sample_qc.ht'
        return sample_qc_path

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import sample_qc

        j = get_batch().new_job(
            'Sample QC',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            if workflow_config.get(key):
                init_batch_args[key] = workflow_config[key]

        j.image(image_path('cpg_workflows'))
        j.command(
            query_command(
                sample_qc,
                sample_qc.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(self.expected_outputs(cohort)),
                str(self.tmp_prefix),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner])
class DenseSubset(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if dense_subset_version := config_retrieve(['large_cohort', 'output_versions', 'dense_subset'], default=None):
            dense_subset_version = slugify(dense_subset_version)

        dense_subset_version = dense_subset_version or get_workflow().output_version
        dense_subset_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / dense_subset_version / 'dense_subset.mt'
        )
        return dense_subset_path

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import dense_subset

        j = get_batch().new_job(
            'Dense Subset',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                dense_subset,
                dense_subset.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, DenseSubset])
class Relatedness(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if relatedness_version := config_retrieve(['large_cohort', 'output_versions', 'relatedness'], default=None):
            relatedness_version = slugify(relatedness_version)

        relatedness_version = relatedness_version or get_workflow().output_version
        relatedness_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version / 'relatedness.ht'
        )
        relatedness_to_drop_path = (
            cohort.analysis_dataset.prefix() / get_workflow().name / relatedness_version / 'relateds_to_drop.ht'
        )

        return dict(
            relatedness=relatedness_path,
            relateds_to_drop=relatedness_to_drop_path,
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job
        from cpg_workflows.large_cohort.relatedness import run

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                out_relatedness_ht_path=self.expected_outputs(cohort)['relatedness'],
                out_relateds_to_drop_ht_path=self.expected_outputs(cohort)['relateds_to_drop'],
                tmp_prefix=self.tmp_prefix,
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, DenseSubset, Relatedness])
class Ancestry(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if ancestry_version := config_retrieve(['large_cohort', 'output_versions', 'ancestry'], default=None):
            ancestry_version = slugify(ancestry_version)

        ancestry_version = ancestry_version or get_workflow().output_version
        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / ancestry_version / 'ancestry'
        return dict(
            scores=prefix / 'scores.ht',
            eigenvalues=prefix / 'eigenvalues.ht',
            loadings=prefix / 'loadings.ht',
            inferred_pop=prefix / 'inferred_pop.ht',
            sample_qc_ht=prefix / 'sample_qc_ht.ht',
        )

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_pca import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                dense_mt_path=inputs.as_path(cohort, DenseSubset),
                sample_qc_ht_path=inputs.as_path(cohort, SampleQC),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort,
                    Relatedness,
                    key='relateds_to_drop',
                ),
                tmp_prefix=self.tmp_prefix,
                out_scores_ht_path=self.expected_outputs(cohort)['scores'],
                out_eigenvalues_ht_path=self.expected_outputs(cohort)['eigenvalues'],
                out_loadings_ht_path=self.expected_outputs(cohort)['loadings'],
                out_inferred_pop_ht_path=self.expected_outputs(cohort)['inferred_pop'],
                out_sample_qc_ht_path=self.expected_outputs(cohort)['sample_qc_ht'],
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[SampleQC, Ancestry, Relatedness])
class AncestryPlots(CohortStage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.out_prefix = get_workflow().web_prefix / 'ancestry'
        self.out_fname_pattern = '{scope}_pc{pci}_{pca_suffix}.{ext}'

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        n_pcs = get_config()['large_cohort']['n_pcs']
        # if there is a pca_plot_name given, add this to the output name
        plot_name = get_config()['large_cohort'].get('pca_plot_name')
        pca_suffix = ''
        if plot_name:
            pca_suffix = plot_name.replace('-', '_')
        return {
            str(pc_num): self.out_prefix
            / self.out_fname_pattern.format(
                scope='dataset',
                pci=pc_num,
                pca_suffix=pca_suffix,
                ext='html',
            )
            for pc_num in range(1, n_pcs)
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.ancestry_plots import run
        from cpg_workflows.large_cohort.dataproc_utils import dataproc_job

        j = dataproc_job(
            job_name=self.__class__.__name__,
            function=run,
            function_path_args=dict(
                out_path_pattern=self.out_prefix / self.out_fname_pattern,
                sample_qc_ht_path=inputs.as_path(cohort, Ancestry, key='sample_qc_ht'),
                scores_ht_path=inputs.as_path(cohort, Ancestry, key='scores'),
                eigenvalues_ht_path=inputs.as_path(cohort, Ancestry, key='eigenvalues'),
                loadings_ht_path=inputs.as_path(cohort, Ancestry, key='loadings'),
                inferred_pop_ht_path=inputs.as_path(
                    cohort,
                    Ancestry,
                    key='inferred_pop',
                ),
                relateds_to_drop_ht_path=inputs.as_path(
                    cohort,
                    Relatedness,
                    key='relateds_to_drop',
                ),
            ),
            depends_on=inputs.get_jobs(cohort),
        )
        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=[Combiner, SampleQC, Relatedness])
class MakeSiteOnlyVcf(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if site_only_version := config_retrieve(['large_cohort', 'output_versions', 'makesiteonly'], default=None):
            site_only_version = slugify(site_only_version)

        site_only_version = site_only_version or get_workflow().output_version

        return {
            'as': cohort.analysis_dataset.prefix() / get_workflow().name / site_only_version / 'as_siteonly.vcf.bgz',
            'as_tbi': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / site_only_version
            / 'as_siteonly.vcf.bgz.tbi',
            'quasi': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / site_only_version
            / 'quasi_siteonly.vcf.bgz',
            'quasi_tbi': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / site_only_version
            / 'quasi_siteonly.vcf.bgz.tbi',
            'ht': cohort.analysis_dataset.prefix() / get_workflow().name / site_only_version / 'siteonly.ht',
            'corrected_mt': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / site_only_version
            / 'filtered_and_corrected.mt',
            'pre_adjusted': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / site_only_version
            / 'siteonly_pre_vcf_adjusted.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import site_only_vcf

        j = get_batch().new_job(
            'MakeSiteOnlyVcf',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )

        sitesvcf_config = config_retrieve('sitesvcf')

        j.image(image_path('cpg_workflows'))
        j.memory(sitesvcf_config.get('memory', '4Gi'))
        j.storage(sitesvcf_config.get('storage', '5Gi'))

        init_batch_args: dict[str, str | int] = {}
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if sitesvcf_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        if 'driver_cores' in sitesvcf_config:
            init_batch_args['driver_cores'] = sitesvcf_config['driver_cores']

        j.command(
            query_command(
                site_only_vcf,
                site_only_vcf.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, SampleQC)),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)['as']),
                str(self.expected_outputs(cohort)['quasi']),
                str(self.expected_outputs(cohort)['ht']),
                str(self.expected_outputs(cohort)['corrected_mt']),
                str(self.expected_outputs(cohort)['pre_adjusted']),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, self.expected_outputs(cohort), [j])


@stage(required_stages=MakeSiteOnlyVcf)
class Vqsr(CohortStage):
    """
    The Vqsr stage performs Variant Quality Score Recalibration (VQSR) and generates a site-only VCF file.
    Additionally, it extracts, edits, and saves the VCF header as a separate file to address parsing issues
    in subsequent stages.

    Why this is necessary:
    -----------------------
    - The VCF file generated by the VQSR process contains an `SB` INFO field with incorrect metadata:
      > ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
      Even though the `SB` field contains lists of integers (e.g., `SB=6,11,2,0`), the metadata incorrectly
      specifies it as a single float. This causes Hail to throw a parsing error:
      > java.lang.NumberFormatException: For input string: "6,11,2,0"

    - To avoid this issue, the header is extracted and modified to correct the `SB` field metadata:
      > ##INFO=<ID=SB,Number=.,Type=Float,Description="Strand Bias">

    - The corrected header is saved as a separate file (`header_siteonly.vqsr.vcf.gz`) so that it can be
      used in the subsequent `LoadVqsr` stage to overwrite the original header when importing the VCF into Hail.

    This ensures that the VCF file can be successfully parsed by Hail without errors, enabling downstream
    processing of the recalibrated variants.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if vqsr_version := config_retrieve(['large_cohort', 'output_versions', 'vqsr'], default=None):
            vqsr_version = slugify(vqsr_version)

        vqsr_version = vqsr_version or get_workflow().output_version
        as_or_quasi = config_retrieve(
            ['large_cohort', 'vqsr_input_vcf'],
            default='quasi',
        )
        return {
            'vcf': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / vqsr_version
            / f'{as_or_quasi}_siteonly.vqsr.vcf.gz',
            'tbi': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / vqsr_version
            / f'{as_or_quasi}_siteonly.vqsr.vcf.gz.tbi',
            'reheadered_header': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / vqsr_version
            / 'header_siteonly.vqsr.txt',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.jobs import vqsr

        vcf_path = inputs.as_path(
            cohort,
            MakeSiteOnlyVcf,
            key=config_retrieve(['large_cohort', 'vqsr_input_vcf'], default='quasi'),
        )
        jobs = vqsr.make_vqsr_jobs(
            b=get_batch(),
            input_siteonly_vcf_path=vcf_path,
            gvcf_count=len(cohort.get_sequencing_groups()),
            out_path=self.expected_outputs(cohort)['vcf'],
            tmp_prefix=self.tmp_prefix,
            use_as_annotations=get_config()['workflow'].get('use_as_vqsr', True),
            intervals_path=get_config()['workflow'].get('intervals_path'),
            job_attrs=self.get_job_attrs(),
        )

        # Adapted from https://github.com/populationgenomics/production-pipelines/blob/e944d7d730be606ab255587a0e02d6c0831361b8/cpg_workflows/jobs/gcnv.py#L364-L373
        outputs = self.expected_outputs(cohort)
        b = get_batch()

        reheader_job = b.new_job(
            'ReheaderVcf',
            (self.get_job_attrs() or {}) | {'tool': 'bcftools'},
        )

        reheader_job.depends_on(*jobs)

        reheader_job.image(image_path('bcftools'))
        reheader_job.storage(config_retrieve(['vqsr_reheader', 'storage'], default='16Gi'))

        vqsr_vcf = b.read_input(outputs['vcf'])

        # pull the header into a temp file
        reheader_job.command(f'bcftools view -h {vqsr_vcf} > {reheader_job.ofile}')

        # sed command to swap Float SB to Integer in-place and allow any length
        reheader_job.command(
            rf"sed -i 's/<ID=SB,Number=1,Type=Float/<ID=SB,Number=.,Type=Float/' {reheader_job.ofile}",
        )

        b.write_output(
            reheader_job.ofile,
            str(outputs['reheadered_header']),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)


@stage(required_stages=[MakeSiteOnlyVcf, Vqsr])
class LoadVqsr(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if load_vqsr_version := config_retrieve(['large_cohort', 'output_versions', 'loadvqsr'], default=None):
            load_vqsr_version = slugify(load_vqsr_version)

        load_vqsr_version = load_vqsr_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / get_workflow().name / load_vqsr_version / 'vqsr.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import load_vqsr

        j = get_batch().new_job(
            'LoadVqsr',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                load_vqsr,
                load_vqsr.run.__name__,
                str(inputs.as_path(cohort, MakeSiteOnlyVcf, key='pre_adjusted')),
                str(inputs.as_path(cohort, Vqsr, key='vcf')),
                str(inputs.as_path(cohort, Vqsr, key='reheadered_header')),
                str(self.expected_outputs(cohort)),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[LoadVqsr])
class ConvertSiteOnlyHTToVcfShards(CohortStage):
    """
    Convert the site-only HT to VCF shards.
    This is used to create the VCF fragments for the VEP stage.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        The output is a directory with the VCF fragments.
        """
        if ht_to_vcf_version := config_retrieve(['large_cohort', 'output_versions', 'ht_to_vcf_version'], default=None):
            ht_to_vcf_version = slugify(ht_to_vcf_version)

        ht_to_vcf_version = ht_to_vcf_version or get_workflow().output_version
        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / ht_to_vcf_version
        return {
            # this will be the write path for fragments of sites-only VCF (header-per-shard)
            'hps_vcf_dir': prefix / 'site_only.vqsr.vcf.bgz',
            # this will be the file which contains the name of all fragments (header-per-shard)
            'hps_shard_manifest': prefix / 'site_only.vqsr.vcf.bgz' / SHARD_MANIFEST,
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import convert_siteonly_ht_to_vcf_shards

        j = get_batch().new_job(
            'ConvertSiteOnlyHTToVcfShards',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(config_retrieve(['workflow', 'driver_image']))

        j.command(
            query_command(
                convert_siteonly_ht_to_vcf_shards,
                convert_siteonly_ht_to_vcf_shards.run.__name__,
                str(inputs.as_path(cohort, LoadVqsr)),
                str(self.expected_outputs(cohort)['hps_vcf_dir'].parent / 'repartitioned'),
                str(self.expected_outputs(cohort)['hps_vcf_dir']),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[ConvertSiteOnlyHTToVcfShards])
class LCAnnotateFragmentedVcfWithVep(CohortStage):
    """
    Annotate VCF with VEP.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        """
        Should this be in tmp? We'll never use it again maybe?
        """
        if vep_version := config_retrieve(['large_cohort', 'output_versions', 'vep'], default=None):
            vep_version = slugify(vep_version)

        vep_version = vep_version or get_workflow().output_version

        return cohort.analysis_dataset.prefix() / get_workflow().name / vep_version / 'vep.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput:
        output = self.expected_outputs(cohort)
        manifest_file: Path = inputs.as_path(
            target=cohort,
            stage=ConvertSiteOnlyHTToVcfShards,
            key='hps_shard_manifest',
        )

        if not manifest_file.exists():
            raise ValueError(
                f'Manifest file {str(manifest_file)} does not exist, '
                f'run the large_cohort workflow with workflows.last_stages=[ConvertSiteOnlyHTToVcfShards]',
            )

        with open(manifest_file, 'r') as f:
            manifest = f.read().strip().splitlines()

        input_vcfs = [to_path(manifest_file.parent / vcf) for vcf in manifest]

        if len(input_vcfs) == 0:
            raise ValueError(f'No VCFs in {manifest_file}')

        vep_jobs = add_vep_jobs(
            b=get_batch(),
            input_vcfs=input_vcfs,
            tmp_prefix=self.tmp_prefix / 'tmp',
            scatter_count=len(input_vcfs),
            out_path=output,
            job_attrs=self.get_job_attrs(),
        )

        return self.make_outputs(cohort, data=output, jobs=vep_jobs)


@stage(required_stages=[Combiner, Relatedness, Ancestry, LoadVqsr, LCAnnotateFragmentedVcfWithVep])
class Frequencies(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if frequencies_version := config_retrieve(['large_cohort', 'output_versions', 'frequencies'], default=None):
            frequencies_version = slugify(frequencies_version)

        frequencies_version = frequencies_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / get_workflow().name / frequencies_version / 'frequencies.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import frequencies

        j = get_batch().new_job(
            'Frequencies',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        if 'driver_cores' in workflow_config:
            init_batch_args['driver_cores'] = workflow_config['driver_cores']

        j.image(image_path('cpg_workflows'))

        j.command(
            query_command(
                frequencies,
                frequencies.run.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, Ancestry, key='sample_qc_ht')),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(inputs.as_path(cohort, LoadVqsr)),
                str(inputs.as_path(cohort, LCAnnotateFragmentedVcfWithVep)),
                str(self.expected_outputs(cohort)),
                init_batch_args=init_batch_args,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage(required_stages=[Combiner, Ancestry, Relatedness])
class GenerateCoverageTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if coverage_version := config_retrieve(['large_cohort', 'output_versions', 'coverage'], default=None):
            coverage_version = slugify(coverage_version)

        coverage_version = coverage_version or get_workflow().output_version
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / coverage_version

        return {
            'group_membership_ht': prefix / f'{sequencing_type}_group_membership.ht',
            'coverage_ht': prefix / f'{sequencing_type}_coverage.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import compute_stats_for_all_sites

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            init_batch_args[batch_key] = 'standard' if not workflow_config.get(config_key) else 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            # Default to 1 core for worker_cores if not specified
            if key == 'worker_cores' and not workflow_config.get('worker_cores'):
                init_batch_args[key] = 1
            else:
                init_batch_args[key] = workflow_config.get(key)

        coverage_table_j = get_batch().new_job(
            'GenerateCoverageTable',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        coverage_table_j.memory(init_batch_args['worker_memory']).cpu(init_batch_args['worker_cores'])
        coverage_table_j.spot(config_retrieve(['workflow', 'preemptible_driver'], False))

        coverage_table_j.image(config_retrieve(['workflow', 'driver_image']))
        coverage_table_j.command(
            query_command(
                compute_stats_for_all_sites,
                compute_stats_for_all_sites.run_coverage.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, Ancestry, key='sample_qc_ht')),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)['group_membership_ht']),
                str(self.expected_outputs(cohort)['coverage_ht']),
                setup_gcp=True,
                init_batch_args=init_batch_args,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[coverage_table_j])


@stage(required_stages=[Combiner, Ancestry, Relatedness])
class GenerateAlleleNumberTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if an_version := config_retrieve(['large_cohort', 'output_versions', 'allele_number'], default=None):
            an_version = slugify(an_version)

        an_version = an_version or get_workflow().output_version
        sequencing_type = config_retrieve(['workflow', 'sequencing_type'])
        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / an_version

        return {
            'group_membership_ht': prefix / f'{sequencing_type}_group_membership.ht',
            'an_ht': prefix / f'{sequencing_type}_allele_number.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import compute_stats_for_all_sites

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            init_batch_args[batch_key] = 'standard' if not workflow_config.get(config_key) else 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            # Default to 1 core for worker_cores if not specified
            if key == 'worker_cores' and not workflow_config.get('worker_cores'):
                init_batch_args[key] = 1
            else:
                init_batch_args[key] = workflow_config.get(key)

        allele_number_j = get_batch().new_job(
            'GenerateAlleleNumber',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        allele_number_j.memory(init_batch_args['worker_memory']).cpu(init_batch_args['worker_cores'])
        allele_number_j.spot(config_retrieve(['workflow', 'preemptible_driver'], False))

        allele_number_j.image(config_retrieve(['workflow', 'driver_image']))
        allele_number_j.command(
            query_command(
                compute_stats_for_all_sites,
                compute_stats_for_all_sites.run_an_calculation.__name__,
                str(inputs.as_path(cohort, Combiner, key='vds')),
                str(inputs.as_path(cohort, Ancestry, key='sample_qc_ht')),
                str(inputs.as_path(cohort, Relatedness, key='relateds_to_drop')),
                str(self.expected_outputs(cohort)['group_membership_ht']),
                str(self.expected_outputs(cohort)['an_ht']),
                setup_gcp=True,
                init_batch_args=init_batch_args,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[allele_number_j])


@stage(required_stages=[LoadVqsr])
class VariantBinnedSummaries(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> Path:
        if var_binned_version := config_retrieve(
            ['large_cohort', 'output_versions', 'var_binned_summaries'],
            default=None,
        ):
            var_binned_version = slugify(var_binned_version)

        var_binned_version = var_binned_version or get_workflow().output_version
        return cohort.analysis_dataset.prefix() / get_workflow().name / var_binned_version / 'binned_summary.ht'

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import variant_binned_summaries

        j = get_batch().new_job(
            'VariantBinnedSummaries',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            if workflow_config.get(config_key):
                init_batch_args[batch_key] = 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            if workflow_config.get(key):
                init_batch_args[key] = workflow_config[key]

        happy_vcf_path = config_retrieve(['large_cohort', 'happy_vcf_path'])
        n_bins = config_retrieve(['large_cohort', 'n_bins'], default=100)  # FIXME default also set in function
        fam_stats_ht_path = config_retrieve(['large_cohort', 'fam_stats_ht_path'], default=None)
        use_truth_sample_concordance = config_retrieve(['large_cohort', 'use_truth_sample_concordance'], default=True)

        j.command(
            query_command(
                variant_binned_summaries,
                variant_binned_summaries.create_binned_summary.__name__,
                str(inputs.as_path(cohort, LoadVqsr)),
                str(happy_vcf_path),
                str(self.expected_outputs(cohort)),
                n_bins,
                fam_stats_ht_path,
                use_truth_sample_concordance,
                setup_gcp=True,
                init_batch_args=init_batch_args,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage()
class JointFrequencyTable(CohortStage):
    """
    Generate a joint frequency table for genome and exome variants.
    This stage outputs the frequency table in Hail Table format.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        if joint_freq_version := config_retrieve(['large_cohort', 'output_versions', 'joint_frequency'], default=None):
            joint_freq_version = slugify(joint_freq_version)
        joint_freq_version = joint_freq_version or get_workflow().output_version

        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / joint_freq_version
        return {
            'joint_freq': prefix / 'joint_frequency.ht',
            'contingency_ht': prefix / 'contingency_table_test.ht',
            'cmh_ht': prefix / 'cmh.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import joint_frequencies

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            init_batch_args[batch_key] = 'standard' if not workflow_config.get(config_key) else 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            # Default to 1 core for worker_cores if not specified
            if key == 'worker_cores' and not workflow_config.get('worker_cores'):
                init_batch_args[key] = 1
            else:
                init_batch_args[key] = workflow_config.get(key)

        j = get_batch().new_job(
            'JointFrequencyTable',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        j.memory(init_batch_args['worker_memory']).cpu(init_batch_args['worker_cores'])

        genome_freq_ht = config_retrieve(['large_cohort', 'output_versions', 'frequencies_genome'], default=None)
        exome_freq_ht = config_retrieve(['large_cohort', 'output_versions', 'frequencies_exome'], default=None)
        genome_all_sites_path = config_retrieve(['large_cohort', 'output_versions', 'genome_all_sites'], default=None)
        exome_all_sites_path = config_retrieve(['large_cohort', 'output_versions', 'exome_all_sites'], default=None)

        j.command(
            query_command(
                joint_frequencies,
                joint_frequencies.run.__name__,
                str(genome_freq_ht),
                str(exome_freq_ht),
                str(genome_all_sites_path),
                str(exome_all_sites_path),
                str(self.expected_outputs(cohort)['joint_freq']),
                str(self.expected_outputs(cohort)['contingency_ht']),
                str(self.expected_outputs(cohort)['cmh_ht']),
                setup_gcp=True,
                init_batch_args=init_batch_args,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage()
class GenerateGeneTable(CohortStage):
    """
    Generate a release-ready gene table of genome and exome variants.
    This stage also outputs an intermediate transcripts table for PrepareBrowserTable stage.
    """

    def expected_outputs(self, cohort: Cohort) -> Path:
        if gene_table_version := config_retrieve(['large_cohort', 'output_versions', 'gene_table'], default=None):
            gene_table_version = slugify(gene_table_version)

        gene_table_version = gene_table_version or get_workflow().output_version

        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / gene_table_version
        return {
            'gene_table': prefix / 'gene_table.ht',
            'transcripts_grch38_base': prefix / 'transcripts_grch38_base.ht',
            'mane_select_transcripts': prefix / 'mane_select_transcripts.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import generate_gene_table

        j = get_batch().new_job(
            'GenerateGeneTable',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        genome_freq_ht = config_retrieve(['large_cohort', 'output_versions', 'frequencies_genome'], default=None)
        exome_freq_ht = config_retrieve(['large_cohort', 'output_versions', 'frequencies_exome'], default=None)

        j.command(
            query_command(
                generate_gene_table,
                generate_gene_table.run.__name__,
                str(genome_freq_ht),
                str(exome_freq_ht),
                str(self.tmp_prefix / 'browser'),
                str(self.expected_outputs(cohort)['transcripts_grch38_base']),
                str(self.expected_outputs(cohort)['mane_select_transcripts']),
                str(self.expected_outputs(cohort)['gene_table']),
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


# @stage(required_stages=[Frequencies])
@stage()
class PrepareBrowserTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        if browser_version := config_retrieve(['large_cohort', 'output_versions', 'preparebrowsertable'], default=None):
            browser_version = slugify(browser_version)

        browser_version = browser_version or get_workflow().output_version
        return {
            'browser': cohort.analysis_dataset.prefix() / get_workflow().name / browser_version / 'browser.ht',
            'exome_variants': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / browser_version
            / 'exome_variants.ht',
            'genome_variants': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / browser_version
            / 'genome_variants.ht',
            'joint_freq': cohort.analysis_dataset.prefix()
            / get_workflow().name
            / browser_version
            / 'joint_freq_release.ht',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort import browser_prepare

        j = get_batch().new_job(
            'PrepareBrowserTable',
            (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
        )
        j.image(image_path('cpg_workflows'))

        exome_freq_ht_path = config_retrieve(['large_cohort', 'output_versions', 'frequencies_exome'], default=None)
        genome_freq_ht_path = config_retrieve(['large_cohort', 'output_versions', 'frequencies_genome'], default=None)
        transcripts_grch38_base_path = config_retrieve(
            ['large_cohort', 'output_versions', 'transcripts_grch38_base'],
            default=None,
        )
        mane_select_transcripts_path = config_retrieve(
            ['large_cohort', 'output_versions', 'mane_select_transcripts'],
            default=None,
        )
        joint_freq_ht_path = config_retrieve(['large_cohort', 'output_versions', 'joint_frequency'], default=None)
        contingency_ht_path = config_retrieve(
            ['large_cohort', 'output_versions', 'contingency_table_test'],
            default=None,
        )
        chm_ht_path = config_retrieve(['large_cohort', 'output_versions', 'cmh'], default=None)

        j.command(
            query_command(
                browser_prepare,
                browser_prepare.prepare_v4_variants.__name__,
                # hard-coding Frequencies tables for now
                exome_freq_ht_path,
                genome_freq_ht_path,
                str(joint_freq_ht_path),
                str(contingency_ht_path),
                str(chm_ht_path),
                str(self.expected_outputs(cohort)['browser']),
                str(self.expected_outputs(cohort)['exome_variants']),
                str(self.expected_outputs(cohort)['genome_variants']),
                str(self.expected_outputs(cohort)['joint_freq']),
                transcripts_grch38_base_path,
                mane_select_transcripts_path,
                setup_gcp=True,
            ),
        )

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=[j])


@stage()
class PrepareBrowserVcfDataDownload(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        data_download_config = config_retrieve(['large_cohort', 'output_versions', 'data_download'], default={})
        browser_vcf_version = data_download_config.get('version', None)
        if browser_vcf_version is None:
            raise ValueError(
                'large_cohort.output_versions.data_download.version must be set in config for PrepareBrowserVcfDataDownload stage',
            )
        browser_vcf_version = slugify(browser_vcf_version)

        prefix = cohort.analysis_dataset.prefix() / get_workflow().name / browser_vcf_version
        data_type: str = data_download_config.get('data_type', 'exomes')
        chroms = data_download_config.get('chroms', [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])
        return {
            **{chrom: prefix / data_type / f'{chrom}_variants.vcf.bgz' for chrom in chroms},
            # **{chrom: prefix / f'{data_type}_{chrom}_variants.vcf.bgz.tbi' for chrom in chroms},
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        from cpg_workflows.large_cohort.scripts import browser_vcf_release

        init_batch_args: dict[str, str | int] = {}
        workflow_config = config_retrieve('workflow')

        # Memory parameters
        for config_key, batch_key in [('highmem_workers', 'worker_memory'), ('highmem_drivers', 'driver_memory')]:
            init_batch_args[batch_key] = 'standard' if not workflow_config.get(config_key) else 'highmem'
        # Cores parameter
        for key in ['driver_cores', 'worker_cores']:
            # Default to 1 core for worker_cores if not specified
            if key == 'worker_cores' and not workflow_config.get('worker_cores'):
                init_batch_args[key] = 1
            else:
                init_batch_args[key] = workflow_config.get(key)

        data_download_config = config_retrieve(['large_cohort', 'output_versions', 'data_download'], default={})
        data_type: str = data_download_config.get('data_type', 'exomes')

        freq_ht_path: str = data_download_config.get('ht_to_export', None)

        exome_freq_ht_path: str = data_download_config.get('frequencies_exome', None)
        genome_freq_ht_path: str = data_download_config.get('frequencies_genome', None)
        vqsr_ht_path: str = data_download_config.get('vqsr_ht', None)

        if data_type == 'joint':
            assert (
                exome_freq_ht_path is not None and genome_freq_ht_path is not None
            ), 'Both exome and genome frequency HT paths must be provided for joint data download.'

        if data_type == 'exomes' or data_type == 'genomes':
            assert vqsr_ht_path is not None, 'VQSR HT path must be provided for exome or genome data download.'

        joint_included: bool = data_download_config.get('joint_included', False)

        jobs = []

        # only schedule a repartition job if this is a genome. For exomes use the original HT path
        if data_type == 'exomes':
            prejob = None
        else:
            ht_path_repartitioned = output_path(f"{data_type}_repartitioned.ht", category="tmp")
            prejob = get_batch().new_job(
                f'RepartitionFrequenciesTable_{data_type}',
                (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
            )
            prejob.image(workflow_config['driver_image'])
            prejob.command(
                query_command(
                    browser_vcf_release,
                    browser_vcf_release.repartition_frequencies_table.__name__,
                    freq_ht_path,
                    data_type,
                    ht_path_repartitioned,
                    setup_gcp=True,
                    init_batch_args=init_batch_args,
                ),
            )
            freq_ht_path = ht_path_repartitioned
            jobs.append(prejob)

        outputs = self.expected_outputs(cohort)
        for chrom, vcf_outpath in outputs.items():
            j = get_batch().new_job(
                f'PrepareBrowserVcfDataDownload_{chrom}_{data_type}',
                (self.get_job_attrs() or {}) | {'tool': HAIL_QUERY},
            )
            j.image(workflow_config['driver_image'])
            j.command(
                query_command(
                    browser_vcf_release,
                    browser_vcf_release.run_browser_vcf_data_download.__name__,
                    freq_ht_path,
                    data_type,
                    chrom,
                    str(vcf_outpath),
                    joint_included,
                    str(vqsr_ht_path) if vqsr_ht_path else None,
                    str(exome_freq_ht_path) if exome_freq_ht_path else None,
                    str(genome_freq_ht_path) if genome_freq_ht_path else None,
                    setup_gcp=True,
                    init_batch_args=init_batch_args,
                ),
            )

            if prejob:
                j.depends_on(prejob)

            jobs.append(j)

        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
