#!/usr/bin/env python3

"""
Batch pipeline to run sample QC.
"""

import logging

import click

from cpg_pipes.images import DRIVER_IMAGE
from cpg_pipes.jobs import fastqc, pedigree
from cpg_pipes.jobs.cram_qc import samtools_stats, verify_bamid, picard_wgs_metrics
from cpg_pipes.pipeline.analysis import CramPath
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.pipeline import stage, Pipeline, PipelineError
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput, DatasetStage

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class FastqcStage(SampleStage):
    """
    Run FastQC on alignment inputs.
    """
    def expected_result(self, sample: Sample):
        """
        Stage is expected to generate a FastQC HTML report, and a zip file for 
        parsing with MuiltiQC.
        """
        folder = sample.dataset.get_bucket(self.pipe) / 'qc' / 'fastqc'
        return {
            'html': folder / (sample.id + '_fastqc.html'),
            'zip': folder / (sample.id + '_fastqc.zip'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using the "fastqc" function implemented in the jobs module
        """
        if not sample.alignment_input:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find read data, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(
                    f'No alignment input found for {sample.id}. '
                    f'Checked: Sequence entry and type=CRAM Analysis entry'
                )

        job = fastqc.fastqc(
            b=self.pipe.b,
            output_html_path=self.expected_result(sample)['html'],
            output_zip_path=self.expected_result(sample)['zip'],
            alignment_input=sample.alignment_input,
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[job]
        )
    

@stage
class CramQC(SampleStage):
    """
    Runs alignment QC.
    """
    def expected_result(self, sample: Sample):
        """
        Expected to generate one QC file per tool, to be parsed with MultiQC.
        MultiQC expectations:
        * Samtools stats file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L652-L654
        * Picard file found by contents, so file name can be any:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L539-L541
        * VerifyBAMID file has to have *.selfSM ending:
          https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L783-L784
        """
        folder = sample.dataset.get_bucket(self.pipe) / 'qc'
        return {
            'samtools_stats': folder / 'samtools_stats' / (sample.id + '.stats'),
            'picard_wgs_metrics': folder / 'picard_wgs_metrics' / (sample.id + '.csv'),
            'verify_bamid': folder / 'verify_bamid' / (sample.id + '.selfSM'),
        }

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        cram_path = sample.analysis_cram_path()
        if not cram_path:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find CRAM analysis, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(f'No CRAM analysis found for {sample.id}')

        samtools_stats_j = samtools_stats(
            b=self.pipe.b,
            cram_path=CramPath(cram_path),
            output_path=self.expected_result(sample)['samtools_stats'],
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        picard_wgs_metrics_j = picard_wgs_metrics(
            b=self.pipe.b,
            cram_path=CramPath(cram_path),
            output_path=self.expected_result(sample)['picard_wgs_metrics'],
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        verify_bamid_j = verify_bamid(
            b=self.pipe.b,
            cram_path=CramPath(cram_path),
            output_path=self.expected_result(sample)['verify_bamid'],
            sample_name=sample.id,
            dataset_name=sample.dataset.name,
        )
        return self.make_outputs(
            sample, 
            data=self.expected_result(sample), 
            jobs=[
                samtools_stats_j,
                verify_bamid_j,
                picard_wgs_metrics_j,
            ]
        )


@stage
class CramSomalierStage(SampleStage):
    """
    Genereate fingerprints from CRAMs for pedigree checks.
    """

    def expected_result(self, sample: Sample):
        """
        Expected to generate the fingerprints file
        """
        return sample.get_cram_path(self.pipe).somalier_path

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput:
        """
        Using a function from the jobs module.
        """
        cram_path = sample.analysis_cram_path()
        if not cram_path:
            if self.pipe.skip_samples_with_missing_input:
                logger.error(f'Could not find CRAM analysis, skipping sample {sample.id}')
                sample.active = False
                return self.make_outputs(sample)  # return empty output
            else:
                raise PipelineError(f'No CRAM analysis found for {sample.id}')

        expected_path = self.expected_result(sample)
        j, _ = pedigree.somalier_extact_job(
            b=self.pipe.b,
            sample=sample,
            gvcf_or_cram_or_bam_path=CramPath(cram_path),
            out_fpath=expected_path,
            overwrite=not self.pipe.check_intermediates,
            depends_on=inputs.get_jobs(),
        )
        return self.make_outputs(sample, data=expected_path, jobs=[j])


@stage(required_stages=CramSomalierStage, forced=True)
class CramPedCheckStage(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints
    """

    def expected_result(self, dataset: Dataset):
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        * (optionally) *.somalier-ancestry.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml#L472-L481
        """
        
        prefix = dataset.get_analysis_bucket(self.pipe) / 'qc' / 'somalier'
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            # 'ancestry': prefix / f'{dataset.name}.somalier-ancestry.tsv',
            'html': dataset.get_web_bucket(self.pipe) / 'qc' / 'somalier.html',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        Checks calls job from the pedigree module
        """
        fp_by_sid = inputs.as_path_by_target(stage=CramSomalierStage)

        j, _, _ = pedigree.add_pedigree_jobs(
            self.pipe.b,
            dataset,
            input_path_by_sid=fp_by_sid,
            overwrite=not self.pipe.check_intermediates,
            out_samples_path=self.expected_result(dataset)['samples'],
            out_pairs_path=self.expected_result(dataset)['pairs'],
            # out_ancestry_path=self.expected_result(dataset)['ancestry'],
            out_html_path=self.expected_result(dataset)['html'],
            web_bucket=dataset.get_web_bucket(self.pipe),
            web_url=self.pipe.web_url,
            tmp_bucket=self.pipe.tmp_bucket,
            depends_on=inputs.get_jobs(),
            dry_run=self.pipe.dry_run,
        )
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


@stage(required_stages=[FastqcStage, CramQC, CramPedCheckStage], forced=True)
class MultiQC(DatasetStage):
    """
    Run MultiQC to summarise all QC.
    """
    def expected_result(self, dataset: Dataset):
        return {
            'html': dataset.get_web_bucket(self.pipe) / 'qc' / 'multiqc.html',
            'json': dataset.get_analysis_bucket(self.pipe) / 'qc' / 'multiqc_data.json'
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        fastqc_zip_d = inputs.as_path_by_target(FastqcStage, id='zip')
        samtools_stats_d = inputs.as_path_by_target(CramQC, id='samtools_stats')
        picard_wgs_metrics_d = inputs.as_path_by_target(CramQC, id='picard_wgs_metrics')
        verify_bamid_d = inputs.as_path_by_target(CramQC, id='verify_bamid')
        somalier_samples = inputs.as_path(dataset, CramPedCheckStage, id='samples')
        somalier_pairs = inputs.as_path(dataset, CramPedCheckStage, id='pairs')

        j = self.pipe.b.new_job('Run MultiQC', {'dataset': dataset.name})
        j.image(DRIVER_IMAGE)
    
        qc_endings = set()
        qc_paths = [somalier_samples, somalier_pairs]
        for sid in dataset.get_sample_ids():
            qc_paths.append(fastqc_zip_d[sid])
            qc_paths.append(samtools_stats_d[sid])
            qc_paths.append(picard_wgs_metrics_d[sid])
            qc_paths.append(verify_bamid_d[sid])
            qc_endings.add(fastqc_zip_d[sid].name.replace(sid, ''))
            qc_endings.add(samtools_stats_d[sid].name.replace(sid, ''))
            qc_endings.add(picard_wgs_metrics_d[sid].name.replace(sid, ''))
            qc_endings.add(verify_bamid_d[sid].name.replace(sid, ''))
    
        file_list_path = dataset.get_tmp_bucket(self.pipe) / 'multiqc-file-list.txt'
        with file_list_path.open('w') as f:
            f.writelines([f'{p}\n' for p in qc_paths])
        file_list = self.pipe.b.read_input(str(file_list_path))

        j.env('GOOGLE_APPLICATION_CREDENTIALS', '/gsa-key/key.json')
        j.command(f'pip install multiqc')
        j.cpu(16)
        j.storage('100G')
        j.command(
            f'gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS'
        )
        j.command(f'mkdir inputs')
        j.command(f'cat {file_list} | gsutil -m cp -I inputs/')
    
        ending_list = ', '.join(f'{ending}' for ending in qc_endings)
        mqc_conf = f'extra_fn_clean_exts: [{ending_list}]'
        j.command(
            f'multiqc inputs -o output -f --fn_as_s_name --cl_config "{mqc_conf}"'
        )
        j.command(f'cp output/multiqc_report.html {j.html}')
        j.command(f'cp output/multiqc_data/multiqc_data.json {j.json}')
    
        self.pipe.b.write_output(j.html, str(self.expected_result(dataset)['html']))
        self.pipe.b.write_output(j.json, str(self.expected_result(dataset)['json']))
        
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


@click.command()
@pipeline_click_options
def main(
    input_datasets: list[str],
    output_version: str,
    **kwargs,
):  # pylint: disable=missing-function-docstring
    assert input_datasets
    title = (
        f'QC: joint call from: {", ".join(input_datasets)}, version {output_version}'
    )

    pipeline = Pipeline(
        name='qc_pipeline',
        description=title,
        input_datasets=input_datasets,
        output_version=output_version,
        **kwargs,
    )
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
