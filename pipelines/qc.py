#!/usr/bin/env python3

"""
Batch pipeline to run WGS QC.
"""

import logging
import click
from cloudpathlib import CloudPath
from hailtop.batch.job import Job

from cpg_pipes.hb.batch import Batch
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.images import DRIVER_IMAGE
from cpg_pipes.jobs import fastqc
from cpg_pipes.jobs.cram_qc import samtools_stats, verify_bamid, picard_wgs_metrics
from cpg_pipes.pipeline.analysis import CramPath
from cpg_pipes.pipeline.cli_opts import pipeline_click_options
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.pipeline import stage, Pipeline, PipelineError
from cpg_pipes.pipeline.sample import Sample
from cpg_pipes.pipeline.stage import SampleStage, StageInput, StageOutput, DatasetStage
from pipelines.somalier import CramSomalierPedigree, CramSomalierAncestry

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@stage
class FastQC(SampleStage):
    """
    Run FastQC on alignment inputs.
    """
    def expected_result(self, sample: Sample):
        """
        Stage is expected to generate a FastQC HTML report, and a zip file for 
        parsing with MuiltiQC.
        """
        folder = sample.dataset.get_bucket() / 'qc'
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
        folder = sample.dataset.get_bucket() / 'qc'
        return {
            'samtools_stats': 
                folder / (sample.id + '_samtools_stats.txt'),
            'picard_wgs_metrics': 
                folder / (sample.id + '_picard_wgs_metrics.csv'),
            'verify_bamid': 
                folder / (sample.id + '_verify_bamid.selfSM'),
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


@stage(
    required_stages=[FastQC, CramQC, CramSomalierPedigree, CramSomalierAncestry], 
    forced=True
)
class MultiQC(DatasetStage):
    """
    Run MultiQC to summarise all QC.
    """
    def expected_result(self, dataset: Dataset):
        return {
            'html': dataset.get_web_bucket() / 'qc' / 'multiqc.html',
            'json': dataset.get_analysis_bucket() / 'qc' / 'multiqc_data.json'
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        fastqc_zip_d = inputs.as_path_by_target(FastQC, id='zip')
        samtools_stats_d = inputs.as_path_by_target(CramQC, id='samtools_stats')
        picard_wgs_metrics_d = inputs.as_path_by_target(CramQC, id='picard_wgs_metrics')
        verify_bamid_d = inputs.as_path_by_target(CramQC, id='verify_bamid')
        somalier_samples = inputs.as_path(dataset, CramSomalierPedigree, id='samples')
        somalier_pairs = inputs.as_path(dataset, CramSomalierPedigree, id='pairs')
        somalier_ancestry = inputs.as_path(dataset, CramSomalierAncestry, id='tsv')

        json_path = self.expected_result(dataset)['json']
        html_path = self.expected_result(dataset)['html']
        html_url = str(html_path).replace(
            str(html_path.parent), dataset.get_web_url()
        )

        paths = [somalier_samples, somalier_pairs, somalier_ancestry]
        ending_to_trim = set()  # endings to trim to get sample names
        for sid in dataset.get_sample_ids():
            paths.append(fastqc_zip_d[sid])
            paths.append(samtools_stats_d[sid])
            paths.append(picard_wgs_metrics_d[sid])
            paths.append(verify_bamid_d[sid])
            ending_to_trim.add(fastqc_zip_d[sid].name.replace(sid, ''))
            ending_to_trim.add(samtools_stats_d[sid].name.replace(sid, ''))
            ending_to_trim.add(picard_wgs_metrics_d[sid].name.replace(sid, ''))
            ending_to_trim.add(verify_bamid_d[sid].name.replace(sid, ''))
        modules_to_trim_endings = {
            'fastqc', 'samtools', 'picard/wgs_metrics', 'verifybamid/selfsm'
        }

        j = multiqc(
            self.pipe.b,
            dataset_name=dataset.name,
            tmp_bucket=dataset.get_tmp_bucket(),
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
        )
        
        return self.make_outputs(dataset, data=self.expected_result(dataset), jobs=[j])


def multiqc(
    b: Batch,
    dataset_name: str,
    tmp_bucket: CloudPath,
    paths: list[CloudPath],
    out_html_path: CloudPath,
    out_json_path: CloudPath,
    out_html_url: str | None = None,
    ending_to_trim: set[str] | None = None,
    modules_to_trim_endings: set[str] | None = None,
) -> Job:
    """
    Run MultiQC for the files in `qc_paths`
    @param b: batch object
    @param dataset_name: dataset name
    @param tmp_bucket: bucket for tmp files
    @param paths: file bucket paths to pass into MultiQC 
    @param out_json_path: where to write MultiQC-generated JSON file
    @param out_html_path: where to write the HTML report
    @param out_html_url: URL corresponding to the HTML report
    @param ending_to_trim: trim these endings from input files to get sample names
    @param modules_to_trim_endings: list of modules for which trim the endings
    @return: 
    """
    j = b.new_job('Run MultiQC', {'dataset': dataset_name})
    j.image(DRIVER_IMAGE)

    file_list_path = tmp_bucket / 'multiqc-file-list.txt'
    with file_list_path.open('w') as f:
        f.writelines([f'{p}\n' for p in paths])
    file_list = b.read_input(str(file_list_path))

    j.env('GOOGLE_APPLICATION_CREDENTIALS', '/gsa-key/key.json')
    j.command(f'pip install multiqc')
    STANDARD.set_resources(j, ncpu=16)
    
    endings_conf = ', '.join(list(ending_to_trim)) if ending_to_trim else ''
    modules_conf = ', '.join(list(modules_to_trim_endings)) if modules_to_trim_endings else ''

    cmd = f"""\
    mkdir inputs
    cat {file_list} | gsutil -m cp -I inputs/

    multiqc inputs -o output -f \\
    --cl_config "extra_fn_clean_exts: [{endings_conf}]" \\
    --cl_config "use_filename_as_sample_name: [{modules_conf}]"
    
    cp output/multiqc_report.html {j.html}
    cp output/multiqc_data/multiqc_data.json {j.json}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'
    j.command(wrap_command(cmd, setup_gcp=True))

    b.write_output(j.html, str(out_html_path))
    b.write_output(j.json, str(out_json_path))
    return j


@click.command()
@pipeline_click_options
def main(
    **kwargs,
):  # pylint: disable=missing-function-docstring

    pipeline = Pipeline(
        name='qc_pipeline',
        description='QC',
        **kwargs,
    )
    pipeline.submit_batch()


if __name__ == '__main__':
    main()  # pylint: disable=E1120
