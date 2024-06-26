"""
Plot ancestry PCA analysis results
"""

import logging
from collections import Counter
from typing import Iterable, List

import numpy as np
import pandas as pd
from bokeh.embed import file_html
from bokeh.models import CategoricalColorMapper, HoverTool
from bokeh.palettes import d3, turbo  # flake: disable=F401
from bokeh.plotting import ColumnDataSource, figure
from bokeh.resources import CDN
from bokeh.transform import factor_cmap, factor_mark

import hail as hl

from cpg_utils import Path
from cpg_utils.config import get_config, reference_path
from cpg_utils.hail_batch import genome_build

PROVIDED_LABEL = 'Provided ancestry'
INFERRED_LABEL = 'Inferred ancestry'


def run(
    out_path_pattern: Path,
    sample_qc_ht_path: Path,
    scores_ht_path: Path,
    eigenvalues_ht_path: Path,
    loadings_ht_path: Path,
    inferred_pop_ht_path: Path,
    relateds_to_drop_ht_path: Path,
):
    """
    Generate plots in HTML format, write for each PC (of n_pcs) and
    scope ("dataset", "superpopulation", "population") plus for loadings into
    file paths defined by `out_path_pattern`.
    """
    sample_ht = hl.read_table(str(sample_qc_ht_path))
    scores_ht = hl.read_table(str(scores_ht_path))
    eigenvalues_ht = hl.read_table(str(eigenvalues_ht_path))
    loadings_ht = hl.read_table(str(loadings_ht_path))
    inferred_pop_ht = hl.read_table(str(inferred_pop_ht_path))
    relateds_to_drop_ht = hl.read_table(str(relateds_to_drop_ht_path))

    def filter_relateds(table, relateds_to_drop):
        return table.filter(hl.is_defined(relateds_to_drop[table.s]), keep=False)

    remove_relateds = get_config()['large_cohort'].get('remove_relateds', False)
    if remove_relateds:
        logging.info('Removing relateds from tables prior to plotting.')
        tables = [sample_ht, scores_ht, inferred_pop_ht]
        sample_ht, scores_ht, inferred_pop_ht = [filter_relateds(table, relateds_to_drop_ht) for table in tables]

    scores_ht = scores_ht.annotate(
        superpopulation=sample_ht[scores_ht.s].superpopulation,
        population=sample_ht[scores_ht.s].population,
        training_pop=inferred_pop_ht[scores_ht.s].training_pop,
        is_training=inferred_pop_ht[scores_ht.s].is_training,
        dataset=sample_ht[scores_ht.s].dataset,
    ).cache()

    def key_by_external_id(ht_, meta_ht=None):
        """
        Assuming ht.s is a CPG id, replaces it with external ID,
        assuming it's defined in meta_ht.external_id. This item is
        configurable in the large_cohort toml, under use_external_id
        where the default is false.
        """
        if meta_ht is None:
            meta_ht = ht_
        ht_ = ht_.annotate(old_s=ht_.s).key_by('old_s')
        ht_ = (
            ht_.annotate(
                s=hl.if_else(
                    hl.is_defined(meta_ht[ht_.old_s]),
                    meta_ht[ht_.old_s].external_id,
                    ht_.old_s,
                ),
            )
            .key_by('s')
            .drop('old_s')
        )
        return ht_

    # Key samples by their external IDs
    use_external_id = get_config()['large_cohort']['use_external_id']
    ht = key_by_external_id(scores_ht, sample_ht) if use_external_id else scores_ht.cache()

    # Use eigenvalues to calculate variance
    eigenvalues = eigenvalues_ht.f0.collect()
    eigenvalues_df = pd.to_numeric(eigenvalues)
    variance = np.divide(eigenvalues_df[1:], float(eigenvalues_df.sum())) * 100
    variance = variance.round(2)
    num_pcs_to_plot = len(eigenvalues_df) - 1

    plots = []

    sample_names = ht.s.collect()
    datasets = ht.dataset.collect()
    use_inferred = get_config()['large_cohort']['pca_background']['inferred_ancestry']
    # if the inferred ancestry is set to true in the config, annotate the PCA with the
    # inferred population ancestry (calculated in the ancestry_pca.py script
    superpopulation_label = ht.training_pop.collect() if use_inferred else ht.superpopulation.collect()
    population_label = ht.training_pop.collect() if use_inferred else ht.population.collect()
    # Change 'none' values to dataset name
    analysis_dataset_name = get_config()['workflow']['dataset']
    # TODO: Input datasets will soon be deprecated, please switch to input_cohorts.
    workflow_dataset = get_config()['workflow'].get('input_datasets', [analysis_dataset_name])
    # join dataset names with underscore, in case there are multiple
    workflow_dataset = '_'.join(workflow_dataset)
    superpopulation_label = [workflow_dataset if x is None else x for x in superpopulation_label]
    population_label = [workflow_dataset if x is None else x for x in population_label]
    is_training = ht.is_training.collect() if use_inferred else [False] * len(ht.is_training.collect())
    for scope, title, labels in [
        ('dataset', 'Dataset', datasets),
        ('superpopulation', 'Superpopulation', superpopulation_label),
        ('population', 'Population', population_label),
    ]:
        plots.extend(
            _plot_pca(
                scope=scope,
                title=title,
                labels=labels,
                number_of_pcs=num_pcs_to_plot,
                variance=variance,
                ht=ht,
                datasets=datasets,
                is_training=is_training,
                sample_names=sample_names,
                out_path_pattern=out_path_pattern,
            ),
        )

    plots.extend(_plot_loadings(num_pcs_to_plot, loadings_ht, out_path_pattern=out_path_pattern))

    return plots


def _plot_pca(
    scope,
    title,
    labels,
    number_of_pcs,
    variance,
    ht,
    datasets,
    sample_names,
    is_training,
    out_path_pattern=None,
):
    cntr = Counter(labels)
    # count the number of samples for each group and add it to the labels
    labels = [f'{x} ({cntr[x]})' for x in labels]
    unique_labels = list(Counter(labels).keys())
    palette = turbo(len(unique_labels))
    # if there is a pca_plot_name given, add this to the output name
    plot_name = get_config()['large_cohort'].get('pca_plot_name')
    pca_suffix = ''
    if plot_name:
        pca_suffix = plot_name.replace('-', '_')

    tooltips = [('labels', '@label'), ('samples', '@samples')]
    plots = []
    for i in range(number_of_pcs - 1):
        pc1 = i
        pc2 = i + 1
        plot = figure(
            title=title,
            x_axis_label=f'PC{pc1 + 1} ({variance[pc1]}%)',
            y_axis_label=f'PC{pc2 + 1} ({variance[pc2]}%)',
            tooltips=tooltips,
            width=1000,
        )
        source = ColumnDataSource(
            dict(
                x=ht.scores[pc1].collect(),
                y=ht.scores[pc2].collect(),
                label=labels,
                samples=sample_names,
                dataset=datasets,
                is_training=[{True: PROVIDED_LABEL, False: INFERRED_LABEL}.get(v) for v in is_training],
            ),
        )
        plot.scatter(
            'x',
            'y',
            alpha=0.5,
            marker=factor_mark('is_training', ['cross', 'circle'], [PROVIDED_LABEL, INFERRED_LABEL]),
            source=source,
            size=5,
            color=factor_cmap('label', palette, unique_labels),
            legend_group='label',
        )
        plot.add_layout(plot.legend[0], 'left')
        plots.append(plot)

        if out_path_pattern:
            html = file_html(plot, CDN, title)
            plot_filename_html = str(out_path_pattern).format(scope=scope, pci=pc2, pca_suffix=pca_suffix, ext='html')
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def _plot_loadings(number_of_pcs, loadings_ht, out_path_pattern=None):
    plots = []
    # if there is a pca_plot_name given, add this to the output name
    plot_name = get_config()['large_cohort'].get('pca_plot_name')
    pca_suffix = ''
    if plot_name:
        pca_suffix = plot_name.replace('-', '_')
    gtf_ht = hl.experimental.import_gtf(
        reference_path('gatk_sv/protein_coding_gtf'),
        reference_genome=genome_build(),
        skip_invalid_contigs=True,
        min_partitions=12,
        force_bgz=True,
    )
    for i in range(number_of_pcs - 1):
        pc = i + 1
        plot = manhattan_loadings(
            iteration=i,
            gtf=gtf_ht,
            loadings=loadings_ht,
            title='Loadings of PC ' + str(pc),
            collect_all=True,
        )
        plots.append(plot)
        if out_path_pattern:
            html = file_html(plot, CDN, 'my plot')
            plot_filename_html = str(out_path_pattern).format(
                scope='loadings',
                pci=pc,
                pca_suffix=pca_suffix,
                ext='html',
            )
            with hl.hadoop_open(plot_filename_html, 'w') as f:
                f.write(html)
    return plots


def manhattan_loadings(
    iteration,
    gtf,
    loadings,
    title=None,
    size=4,
    hover_fields=None,
    collect_all=False,
    n_divisions=500,
):
    """modify hail manhattan plot"""
    palette = [
        '#1f77b4',
        '#ff7f0e',
        '#2ca02c',
        '#d62728',
        '#9467bd',
        '#8c564b',
        '#e377c2',
        '#7f7f7f',
        '#bcbd22',
        '#17becf',
    ]

    # add gene names, p-values, and locus info
    loadings = loadings.annotate(gene_names=gtf[loadings.locus].gene_name)
    pvals = hl.abs(loadings.loadings[iteration])
    locus = loadings.locus

    if hover_fields is None:
        hover_fields = {}

    hover_fields['locus'] = hl.str(locus)
    hover_fields['gene'] = hl.str(loadings.gene_names)

    source_pd = hl.plot.plots._collect_scatter_plot_data(  # pylint: disable=protected-access
        ('_global_locus', locus.global_position()),
        ('_pval', pvals),
        fields=hover_fields,
        n_divisions=None if collect_all else n_divisions,
    )
    source_pd['p_value'] = source_pd['_pval']
    source_pd['_contig'] = [locus.split(':')[0] for locus in source_pd['locus']]

    observed_contigs = set(source_pd['_contig'])
    ref = locus.dtype.reference_genome
    observed_contigs = [contig for contig in ref.contigs.copy() if contig in observed_contigs]

    contig_ticks = [
        ref._contig_global_position(contig) + ref.contig_length(contig) // 2  # pylint: disable=protected-access
        for contig in observed_contigs
    ]
    color_mapper = CategoricalColorMapper(factors=ref.contigs, palette=palette[:2] * int((len(ref.contigs) + 1) / 2))

    p = figure(title=title, x_axis_label='Chromosome', y_axis_label='Loadings', width=1000)
    (
        p,
        _,
        legend,
        _,
        _,
        _,
    ) = hl.plot.plots._get_scatter_plot_elements(  # pylint: disable=protected-access
        p,
        source_pd,
        x_col='_global_locus',
        y_col='_pval',
        label_cols=['_contig'],
        colors={'_contig': color_mapper},
        size=size,
    )
    legend.visible = False
    p.xaxis.ticker = contig_ticks
    p.xaxis.major_label_overrides = dict(zip(contig_ticks, observed_contigs))
    p.select_one(HoverTool).tooltips = [t for t in p.select_one(HoverTool).tooltips if not t[0].startswith('_')]

    return p
