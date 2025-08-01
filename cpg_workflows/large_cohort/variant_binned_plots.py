import logging
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd

import hail as hl
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, reference_path
from cpg_workflows.utils import can_reuse, exists
from gnomad.variant_qc.evaluation import compute_grouped_binned_ht
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg


def plot_binned_metric(
    df,
    bin_id,
    y_axis_label,
    y_func,
    x_axis_label="Bin",
    x_func="bin",
    variant_type="SNPs",
    bin_threshold=None,
    tooltips=['n', 'min_score', 'max_score'],
    cumulative=False,
):
    """
    Plot a specific variant metric against QC score ranking bins, to enable the selection of a bin sensitivity threshold.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with 'bin_id', 'bin', 'snv' and metric columns
    bin_id : str
        One of: 'adj_biallelic_bin', 'adj_biallelic_singleton_bin', 'adj_bin',
        'adj_singleton_bin', 'biallelic_bin', 'biallelic_singleton_bin', 'bin', 'singleton_bin'
    y_axis_label : str
        Label for the y-axis
    y_func : str or callable
        Either a column name (e.g., 'n') or a function that takes the DataFrame
        and returns a Series (e.g., lambda x: x['n_ti'] / x['n_tv'])
    x_axis_label : str, optional
        Label for the x-axis. Default is 'Bin'.
    x_func : str or callable, optional
        Either a column name (e.g., 'n') or a function that takes the DataFrame. Default is to use the bin.
    variant_type : str, optional
        One of 'SNPs' or 'Indels', indicating the type of variant data to plot. Default is 'SNPs'.
    bin_threshold : float, optional
        If provided, draws a vertical red dotted line at this bin value
    tooltips : list, optional
        List of column names to include in hover tooltips. Default is ['n','min_score','max_score']
    cumulative : bool, optional
        If True, plot cumulative values from bin 1 to each bin. Default is False.

    Returns:
    --------
    p : bokeh.plotting.figure
        Bokeh figure object
    """

    if variant_type not in ["SNPs", "Indels"]:
        raise ValueError(f"Variant type must be either 'SNPs' or 'Indels', got '{variant_type}'")
    snv = variant_type == "SNPs"  # True if variant_type == SNPs

    if cumulative and x_func != "bin":
        raise ValueError("Cumulative plots can only be generated when 'bin' is on the x-axis")

    # Filter and prepare data
    filtered_df = df[(df['bin_id'] == bin_id) & (df['snv'] == snv)].copy()

    if filtered_df.empty:
        raise ValueError(f"No data points found for bin_id='{bin_id}' and variant_type={variant_type}")

    filtered_df = filtered_df.sort_values('bin')

    # Apply cumulative transformation if requested
    if cumulative:
        numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns.drop('bin')
        filtered_df[numeric_cols] = filtered_df[numeric_cols].cumsum()

    # Build title
    subtitle = f"({bin_id.replace('bin', '').replace('_', ' ').title()}{variant_type})".replace("( ", "(")
    title = f"{y_axis_label} {subtitle}"
    if cumulative:
        title += ' - cumulative'

    # Get axis values
    x_values = filtered_df[x_func] if isinstance(x_func, str) else x_func(filtered_df)
    y_values = filtered_df[y_func] if isinstance(y_func, str) else y_func(filtered_df)

    # Prepare data source
    source_data = {
        'x': x_values,
        'y': y_values,
        **{col: filtered_df[col] for col in tooltips if col in filtered_df.columns},
    }
    source = ColumnDataSource(data=source_data)

    # Create plot
    p = figure(
        width=800,
        height=500,
        title=title,
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        # x_range=(0, 101)
    )

    # Add scatter plot
    scatter = p.scatter(
        x='x',
        y='y',
        source=source,
        size=6,
        alpha=0.8,
        color='navy',
        hover_fill_color='firebrick',
        hover_alpha=1.0,
        hover_line_color='darkred',
    )

    # Configure hover tool
    hover_tooltips = [(x_axis_label, '@x{0.00}'), (y_axis_label, '@y{0.00}')] + [
        (col, f'@{col}') for col in tooltips if col in source_data
    ]

    p.add_tools(HoverTool(renderers=[scatter], tooltips=hover_tooltips))

    # Add threshold line if specified
    if bin_threshold is not None:
        p.add_layout(
            Span(
                location=bin_threshold,
                dimension='height',
                line_color='red',
                line_dash='dashed',
                line_width=2,
                line_alpha=0.7,
            ),
        )

    # Style the plot
    p.title.text_font_size = '14pt'
    p.title.text_font_style = 'bold'
    p.xaxis.axis_label_text_font_size = '12pt'
    p.yaxis.axis_label_text_font_size = '12pt'
    p.grid.grid_line_alpha = 0.3

    return p


def plot_metric_tabs(
    df,
    bin_ids,
    variant_types,
    y_axis_label,
    y_func,
    x_axis_label="Bin",
    x_func="bin",
    adjusted=False,
    snp_bin_threshold=None,
    indel_bin_threshold=None,
    cumulative=[False, True],
    tooltips=['n', 'min_score', 'max_score'],
    width=600,
    height=400,
):
    """
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with 'bin_id', 'bin', 'snv' and metric columns
    bin_ids : list of str
        List of bin_id values to create tabs for
    variant_types : list of variant types to plot
        Allowed variant types are 'SNPs' and/or 'Indels'
    y_axis_label : str
        Label for the y-axis
    y_func : str or callable
        Either a column name or a function that returns values
    x_axis_label : str, optional
        Label for the x-axis. Default is 'Bin'.
    x_func : str or callable, optional
        Either a column name (e.g., 'n') or a function that takes the DataFrame. Default is to use the bin.
    adjusted : bool, optional
        If True, prepends 'adj_' to each of the bin_ids (default: False)
    bin_threshold : float, optional
        If provided, draws a vertical red dotted line at this bin value
    tooltips : list, optional
        List of column names to include in hover tooltips
    width : int, optional
        Plot width in pixels (default: 600)
    height : int, optional
        Plot height in pixels (default: 400)

    Returns:
    --------
    tabs : bokeh.models.Tabs
        Bokeh Tabs object containing all the tabbed plots
    """

    tabs = []

    # Select adjusted genotype bins if requested
    if adjusted:
        bin_ids = [f'adj_{bin_id}' for bin_id in bin_ids]

    for bin_id in bin_ids:
        plots = []
        for variant_type in variant_types:
            try:
                row_plots = []
                for cumulative_plot in cumulative:
                    p = plot_binned_metric(
                        df=df,
                        bin_id=bin_id,
                        y_axis_label=y_axis_label,
                        y_func=y_func,
                        x_axis_label=x_axis_label,
                        x_func=x_func,
                        variant_type=variant_type,
                        bin_threshold=snp_bin_threshold if variant_type == "SNPs" else indel_bin_threshold,
                        tooltips=tooltips,
                        cumulative=cumulative_plot,
                    )
                    p.width = width
                    p.height = height
                    row_plots.append(p)

                plots.append(row_plots)

            except ValueError as e:
                print(f"Warning: Skipping bin_id '{bin_id}' for variant type={variant_type} - {str(e)}")

        if plots:
            if len(plots) == 1:
                layout = row(*plots[0])
            else:
                layout = gridplot(plots, toolbar_location='above')

            tab_title = bin_id.replace('_', ' ').title().replace(' Bin', '').replace('Bin', 'All')
            tab_panel = TabPanel(child=layout, title=tab_title)
            tabs.append(tab_panel)

    if not tabs:
        raise ValueError("No valid plots could be created for any of the provided bin_ids")

    return Tabs(tabs=tabs)


def run(
    binned_summary_ht_path: str,
    binned_plots_outpath: str,
    snp_bin_threshold: int,
    indel_bin_threshold: int,
) -> None:
    """
    Output plots from the binned_summaries_ht table

    :param binned_summary_ht_path: Input Hail table of variant binned summaries
    :param binned_plots_outpath: Output directory for plots
    :param snp_bin_threshold: QC threshold for SNPs
    :param indel_bin_threshold: QC threshold for INDELs
    :return: Outputs plots sequentially
    """
    # Read in the binned_summary_ht table and convert to pandas.
    binned_summary_ht = hl.read_table(binned_summary_ht_path)
    df = binned_summary_ht.to_pandas()  # FIXME refactor df to appropriate name

    # Ti/Tv ratio
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Ti/Tv',
        y_func=lambda df: df['n_ti'] / df['n_tv'],
        snp_bin_threshold=snp_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_ti', 'n_tv'],
    ).save(f"{binned_plots_outpath}/ti_tv.html")

    # Proportion singletons
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin'],
        adjusted=False,
        y_axis_label='Proportion singletons',
        y_func=lambda df: df['n_singleton'] / df['n'],
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_singleton'],
    ).save(f"{binned_plots_outpath}/proportion_singletons.html")

    # Proportion singletons adjusted
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin'],
        adjusted=True,
        y_axis_label='Proportion singletons',
        y_func=lambda df: df['n_singleton'] / df['n'],
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_singleton'],
    ).save(f"{binned_plots_outpath}/proportion_singletons_adj.html")

    # Proportion bi-allelic
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'singleton_bin'],
        adjusted=False,
        y_axis_label='Proportion bi-allelic',
        y_func=lambda df: df['n_biallelic'] / df['n'],
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_biallelic'],
    ).save(f"{binned_plots_outpath}/biallelics.html")

    # ClinVar variants
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Clinvar',
        y_func='n_clinvar',
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_clinvar'],
    ).save(f"{binned_plots_outpath}/clinvar.html")

    # ClinVar pathogenic variants
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Clinvar Pathogenic / Likely Pathogenic',
        y_func='n_clinvar_path',
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_clinvar_path'],
    ).save(f"{binned_plots_outpath}/clinvar_path.html")

    # INDEL ratios
    plot_metric_tabs(
        df=df,
        variant_types=['Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Ins/Del',
        y_func=lambda df: df['n_ins'] / df['n_del'],
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_ins', 'n_del'],
    ).save(f"{binned_plots_outpath}/indel_ratios.html")

    # INDEL 1bp ratios
    plot_metric_tabs(
        df=df,
        variant_types=['Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='1bp Ins/Del',
        y_func=lambda df: df['n_1bp_ins'] / df['n_1bp_del'],
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_1bp_ins', 'n_1bp_del'],
    ).save(f"{binned_plots_outpath}/indel_1bp_ratios.html")

    # INDEL 2bp ratios
    plot_metric_tabs(
        df=df,
        variant_types=['Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='2bp Ins/Del',
        y_func=lambda df: df['n_2bp_ins'] / df['n_2bp_del'],
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_2bp_ins', 'n_2bp_del'],
    ).save(f"{binned_plots_outpath}/indel_2bp_ratios.html")

    # INDEL 3bp ratios
    plot_metric_tabs(
        df=df,
        variant_types=['Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='3bp Ins/Del',
        y_func=lambda df: df['n_3bp_ins'] / df['n_3bp_del'],
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_3bp_ins', 'n_3bp_del'],
    ).save(f"{binned_plots_outpath}/indel_3bp_ratios.html")

    # Truth sample precision
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Precision',
        y_func=lambda df: df['n_tp'] / (df['n_tp'] + df['n_fp']),
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_tp', 'n_fp', 'n_fn'],
        cumulative=[False],
    ).save(f"{binned_plots_outpath}/truth_sample_precision.html")

    # Truth sample recall
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Recall',
        y_func=lambda df: df['n_tp'] / (df['n_tp'] + df['n_fn']),
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_tp', 'n_fp', 'n_fn'],
        cumulative=[False],
    ).save(f"{binned_plots_outpath}/truth_sample_recall.html")

    # Truth sample precision x recall
    plot_metric_tabs(
        df=df,
        variant_types=['SNPs', 'Indels'],
        bin_ids=['bin', 'biallelic_bin', 'singleton_bin', 'biallelic_singleton_bin'],
        adjusted=False,
        y_axis_label='Precision',
        y_func=lambda df: df['n_tp'] / (df['n_tp'] + df['n_fp']),
        x_axis_label='Recall',
        x_func=lambda df: df['n_tp'] / (df['n_tp'] + df['n_fn']),
        snp_bin_threshold=snp_bin_threshold,
        indel_bin_threshold=indel_bin_threshold,
        tooltips=['n', 'min_score', 'max_score', 'n_tp', 'n_fp', 'n_fn', 'bin'],
        cumulative=[False],
    ).save(f"{binned_plots_outpath}/truth_sample_precision_x_recall.html")
