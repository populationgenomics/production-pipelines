"""
Driver for computing structural variants from GATK-SV.
Uses metadata from the sample-metadata server.

Example running:

```
pipelines/gatk_sv.py 
    -n test \
    --analysis-dataset acute-care
    --dataset acute-care
```

OR with a config:

```
pipelines/gatk_sv.py --config pipelines/configs/acute-care-test.yml
```

WHERE pipelines/configs/acute-care-test.yml:

```
namespace: test
analysis_dataset: acute-care
datasets: [acute-care]
```
"""

import logging
from typing import Any

import click

from cpg_pipes import Path
from cpg_pipes.pipeline import stage, SampleStage, StageOutput, DatasetStage, \
    CohortStage, pipeline_click_options, create_pipeline
from cpg_pipes.pipeline.pipeline import StageInput
from cpg_pipes.targets import Sample, Dataset, Cohort

from .utils import add_gatksv_job, get_references, get_dockers, SV_CALLERS

logger = logging.getLogger(__file__)


@stage
class GatherSampleEvidence(SampleStage):
    """
    https://github.com/broadinstitute/gatk-sv#gathersampleevidence    
    """
    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Expectes to produce coverage counts, and VCF for each variant caller
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'coverage_counts': 'coverage_counts.tsv.gz'
        }
        for caller in SV_CALLERS:
            fname_by_key[f'{caller}_vcf'] = f'{caller}.vcf.gz'
            fname_by_key[f'{caller}_index'] = f'{caller}.vcf.gz.tbi'

        # Split reads (SR) and Discordant read pairs (PE)
        for k in ['disc', 'split']:
            fname_by_key[f'pesr_{k}'] = f'pesr_{k}.txt.gz'
            fname_by_key[f'pesr_{k}_index'] = f'pesr_{k}.txt.gz.tbi'

        for key, fname in fname_by_key.items():
            stage_name = self.name.lower()
            d[key] = (
                sample.dataset.get_bucket() / 'gatk_sv' / stage_name /
                (sample.id + '-' + fname)
            )
        return d

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        cram_path = sample.dataset.get_bucket() / 'cram' / f'{sample.id}.cram'

        input_dict: dict[str, Any] = {
            'bam_or_cram_file': str(cram_path),
            'bam_or_cram_index': str(cram_path) + '.crai',
            'sample_id': sample.id,
            # 'revise_base_cram_to_bam': True,
        }

        input_dict.update(get_dockers([
            'sv_pipeline_docker',
            'sv_base_mini_docker',
            'samtools_cloud_docker',
            'gatk_docker',
            'genomes_in_the_cloud_docker',
            'cloud_sdk_docker',
            'wham_docker',
            'manta_docker',
            'sv_pipeline_base_docker',
            'gatk_docker_pesr_override',
        ]))

        input_dict.update(get_references([
            'primary_contigs_fai',
            'primary_contigs_list',
            'reference_fasta',
            'reference_index',
            'reference_dict',
            'reference_version',
            'preprocessed_intervals',
            'manta_region_bed',
            'wham_include_list_bed_file',
        ]))

        expected_d = self.expected_outputs(sample)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset=sample.dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
            sample_id=sample.id,
        )

        return self.make_outputs(sample, data=output_dict, jobs=[j])


@stage(required_stages=GatherSampleEvidence)
class EvidenceQC(DatasetStage):
    """
    # https://github.com/broadinstitute/gatk-sv#evidenceqc
    """
    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Expected to return a bunch of batch-level summary files.
        """
        d: dict[str, Path] = dict()
        fname_by_key = {
            'ploidy_matrix':       'ploidy_matrix.bed.gz',
            'ploidy_plots':        'ploidy_plots.tar.gz',
            'WGD_dist':            'WGD_score_distributions.pdf',
            'WGD_matrix':          'WGD_scoring_matrix_output.bed.gz',
            'WGD_scores':          'WGD_scores.txt.gz',
            'bincov_median':       'RD.txt.gz',
            'bincov_matrix':       'RD.txt.gz',
            'bincov_matrix_index': 'RD.txt.gz.tbi',
        }
        for caller in SV_CALLERS:
            for k in ['low', 'high']:
                fname_by_key[f'{caller}_qc_{k}'] = f'{caller}_QC.outlier.{k}'

        for key, fname in fname_by_key.items():
            d[key] = dataset.get_bucket() / 'gatk_sv' / self.name.lower() / fname
        return d

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Queue jobs
        """
        d = inputs.as_dict_by_target(GatherSampleEvidence)

        sids = [s.id for s in dataset.get_samples()]

        input_dict: dict[str, Any] = {
            'batch': dataset.name,
            'samples': sids,
            'run_vcf_qc': True,  # generates <caller>_qc_low/<caller>_qc_high
            'counts': [str(d[sid]['coverage_counts']) for sid in sids],
        }
        for caller in SV_CALLERS:
            input_dict[f'{caller}_vcfs'] = [
                str(d[sid][f'{caller}_vcf']) for sid in sids
            ]

        input_dict |= get_dockers([
            'sv_base_mini_docker',
            'sv_base_docker',
            'sv_pipeline_docker',
            'sv_pipeline_qc_docker',
        ])
        
        input_dict |= get_references([
            'genome_file',
            'wgd_scoring_mask',
        ])

        expected_d = self.expected_outputs(dataset)

        output_dict, j = add_gatksv_job(
            batch=self.b,
            dataset=dataset,
            wfl_name=self.name,
            input_dict=input_dict,
            expected_out_dict=expected_d,
        )

        return self.make_outputs(dataset, data=output_dict, jobs=[j])


@click.command()
@pipeline_click_options
def main(
    **kwargs,
): 
    """
    GATK-SV workflow.
    """
    pipeline = create_pipeline(
        name='gatk_sv',
        **kwargs
    )
    pipeline.run()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
