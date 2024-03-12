#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script outputs a scatterplot plotting chrX and chrY ploidy, and labelling the points by a column in the
Sample QC Table (eg: inferred karyotypic sex, or dataset) from the SampleQC Hail Table.

 analysis-runner --dataset "bioheart" \
    --description "sex ploidy plotter" \
    --access-level "test" \
    --output-dir "hoptan-str/sex_ploidy_plot/dataset" \
    sex_ploidy_plotter.py --file-path=gs://cpg-bioheart-test/large_cohort/1-0/sample_qc.ht/ \
    --label=dataset

    Script taken from: https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/sex_ploidy/sex_ploidy_plotter.py
"""

import click
import hail as hl
from bokeh.plotting import output_file, save
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def sex_ploidy_plotter(file_path, gcs_path, label):
    init_batch()
    sample_qc_table = hl.read_table(file_path)

    # Create a scatterplot matrix using Hail's plotting functions
    p = hl.plot.scatter(
        x=sample_qc_table.chrX_ploidy,
        y=sample_qc_table.chrY_ploidy,
        label=getattr(sample_qc_table, label),
        title='Scatterplot of chrY_ploidy vs chrX_ploidy',
        xlabel='chrX_ploidy',
        ylabel='chrY_ploidy',
    )

    # Save the plot to a local file, then hadoop_copy to copy to GCS bucket
    output_file('local_plot.html')
    save(p)
    hl.hadoop_copy('local_plot.html', gcs_path)


@click.option(
    '--file-path',
    help='GCS file path to SampleQC Hail Table',
    type=str,
)
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='20G'
)
@click.option('--job-memory', help='Memory of the Hail batch job eg 64G', default='8G')
@click.option(
    '--label', help='Column field in SampleQC table eg dataset or sex_karyotype'
)
@click.option(
    '--out-prefix',
    help='Prefix for output file',
)
@click.command()
def main(file_path, job_storage, job_memory, label, out_prefix):
    # Initialise batch
    b = get_batch()

    # Initialise job
    j = b.new_python_job(name=f'Sex ploidy plotter job')
    j.memory(job_memory)
    j.storage(job_storage)

    gcs_output_path = output_path(f'{out_prefix}_sex_ploidy_plot.html', 'analysis')
    j.call(sex_ploidy_plotter, file_path, gcs_output_path, label)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
