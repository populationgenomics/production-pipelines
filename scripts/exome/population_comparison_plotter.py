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
from bokeh.plotting import output_file, save

import hail as hl

from cpg_utils.hail_batch import get_batch, init_batch, output_path


def pop_comparison(sample_qc_path, inferred_pop_path, gcs_path, label):
    init_batch()
    inferred_pop = hl.read_table(inferred_pop_path)
    sample_qc_table = hl.read_table(sample_qc_path)

    # Add a new column to the table that indicates whether the population was correctly inferred
    sample_qc_table = sample_qc_table.annotate(inferred=inferred_pop[sample_qc_table.s].pop)

    # Filter the table to include only the incorrectly inferred populations
    incorrect_inferences = sample_qc_table.filter(~sample_qc_table.inferred)

    # Create a scatterplot matrix using Hail's plotting functions
    p = hl.plot.scatter(
        x=incorrect_inferences.superpopulation,
        y=inferred_pop.pop,
        label=getattr(incorrect_inferences, label),
        title='Scatterplot of Incorrectly Inferred Populations',
        xlabel='Labeled population',
        ylabel='Inferred population',
    )

    # Save the plot to a local file, then hadoop_copy to copy to GCS bucket
    output_file('local_plot.html')
    save(p)
    hl.hadoop_copy('local_plot.html', gcs_path)


@click.option(
    '--sample-qc-path',
    help='GCS file path to SampleQC Hail Table',
    type=str,
)
@click.option(
    '--inferred-pop-path',
    help='GCS file path to Inferred Population Hail Table',
    type=str,
)
@click.option(
    '--job-storage',
    help='Storage of the Hail batch job eg 30G',
    default='20G',
)
@click.option('--job-memory', help='Memory of the Hail batch job eg 64G', default='8G')
@click.option(
    '--label',
    help='Column field in SampleQC table eg dataset or sex_karyotype',
)
@click.option(
    '--out-prefix',
    help='Prefix for output file',
)
@click.command()
def main(sample_qc_path, inferred_pop_path, job_storage, job_memory, label, out_prefix):
    # Initialise batch
    b = get_batch()

    # Initialise job
    j = b.new_python_job(name='Sex ploidy plotter job')
    j.memory(job_memory)
    j.storage(job_storage)

    gcs_output_path = output_path(f'{out_prefix}_sex_ploidy_plot.html', 'analysis')
    j.call(pop_comparison, sample_qc_path, gcs_output_path, label)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
