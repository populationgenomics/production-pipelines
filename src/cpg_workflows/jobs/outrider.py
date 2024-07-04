"""
Perform outlier gene expression analysis with Outrider.
"""

from os.path import basename
from textwrap import dedent

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse


class Outrider:
    """
    Construct an outrider command for performing outlier gene expression analysis.
    """

    def __init__(
        self,
        input_counts: list[str | Path],
        output: hb.ResourceGroup,
        gtf_file: str | Path | None = None,
        nthreads: int = 8,
        pval_cutoff: float = 0.05,
        z_cutoff: float = 0.0,
    ) -> None:
        self.input_counts = input_counts
        assert isinstance(
            self.input_counts,
            list,
        ), f'input_counts must be a list, instead got {type(self.input_counts)}: {self.input_counts}'
        self.input_counts_r_str = ', '.join([f'"{str(f)}"' for f in self.input_counts])
        self.gtf_file_path = str(gtf_file)
        self.output = output
        self.nthreads = nthreads
        self.pval_cutoff = pval_cutoff
        self.z_cutoff = z_cutoff

        # Build OUTRIDER command
        self.command = """
        R --vanilla <<EOF
        library(OUTRIDER)
        library(tidyverse)

        # Create directories
        dir.create("plots", showWarnings = FALSE)
        dir.create("plots/heatmaps", showWarnings = FALSE)
        dir.create("plots/sig_genes", showWarnings = FALSE)
        dir.create("plots/volcano", showWarnings = FALSE)
        dir.create("plots/stats", showWarnings = FALSE)
        dir.create("results", showWarnings = FALSE)
        """
        self.command += f"""
        # Set significance values
        pval_cutoff <- {str(self.pval_cutoff)}
        z_cutoff <- {str(self.z_cutoff)}
        n_parallel_workers <- {str(self.nthreads - 1)}

        input_counts_files <- c({self.input_counts_r_str})
        """
        self.command += """
        input_counts_tbl_list <- lapply(input_counts_files, function(x) {
            read_tsv(x, col_names = TRUE, comment = "#") %>%
                rename_with(~ gsub("\\\\\\\\.bam$", "", basename(.x)), .cols = last_col()) %>%
                select(c(Geneid, last_col()))
        })
        input_counts_tbl <- reduce(input_counts_tbl_list, left_join, by = "Geneid")
        input_counts_df <- input_counts_tbl %>%
            column_to_rownames("Geneid")

        # Create OutriderDataSet object
        ods <- OutriderDataSet(countData = input_counts_df)

        # Save original data
        ods_orig <- ods

        # Read in gene annotations
        """
        self.command += f"""
        gtf_file <- "{self.gtf_file_path}"
        """
        self.command += """
        ods <- filterExpression(
            ods,
            gtfFile = gtf_file,
            filterGenes = FALSE,
            savefpkm = TRUE
        )

        # Plotting
        png(file = "plots/stats/fpkm.png", width = 4000, height = 4000, res = 600)
        plotFPKM(ods)
        dev.off()

        png(file = "plots/stats/expressed_genes_per_sample.png", width = 4000, height = 4000, res = 600)
        plotExpressedGenes(ods)
        dev.off()

        # Save filtered data
        ods_filt <- ods

        # Subset based on filter
        ods <- ods[mcols(ods)\\$passedFilter, ]

        # Plotting - heatmaps
        # also annotates the clusters resulting from the dendrogram
        png(file = "plots/heatmaps/sample_cor_heatmap.png", width = 4000, height = 4000, res = 600)
        ods <- plotCountCorHeatmap(
            ods,
            colGroups = NA,
            rowGroups = NA,
            normalized = FALSE,
            nRowCluster = NA,
            nColCluster = NA
        )
        dev.off()

        # heatmap of the gene/sample expression
        png(file = "plots/heatmaps/gene_sample_heatmap.png", width = 4000, height = 4000, res = 600)
        ods <- plotCountGeneSampleHeatmap(
            ods,
            colGroups = NA,
            rowGroups = NA,
            normalized = FALSE,
            nRowCluster = NA,
            nColCluster = NA
        )
        dev.off()

        # Save subsetted data
        ods_subset <- ods

        # Estimate size factors
        ods <- estimateSizeFactors(ods)

        # Set sample exclusion mask (if necessary)
        sampleExclusionMask(ods) <- FALSE

        # Find optimal encoding dimension
        ods <- findEncodingDim(ods, BPPARAM=MulticoreParam(workers = n_parallel_workers))
        png(file = "plots/stats/enc_dim_search.png", width = 4000, height = 4000, res = 600)
        p <- plotEncDimSearch(ods)
        print(p)
        dev.off()
        optimal_q <- ods@metadata\\$optimalEncDim

        # Save dataset with optimal encoding dimension
        ods_opt_q <- ods

        # Use the autoencoder method to control for confounders
        n_iter <- 15
        ods <- controlForConfounders(ods, q = optimal_q, iterations = n_iter, BPPARAM=MulticoreParam(workers = n_parallel_workers))

        # Plotting - normalized heatmaps
        png(file = "plots/heatmaps/sample_cor_heatmap_normalized.png", width = 4000, height = 4000, res = 600)
        ods <- plotCountCorHeatmap(
            ods,
            colGroups = NA,
            rowGroups = NA,
            normalized = TRUE,
            nRowCluster = NA,
            nColCluster = NA
        )
        dev.off()

        # heatmap of the gene/sample expression
        png(file = "plots/heatmaps/gene_sample_heatmap_normalized.png", width = 4000, height = 4000, res = 600)
        ods <- plotCountGeneSampleHeatmap(
            ods,
            colGroups = NA,
            rowGroups = NA,
            normalized = TRUE,
            nRowCluster = NA,
            nColCluster = NA
        )
        dev.off()

        # Save dataset after controlling for confounders
        ods_controlled <- ods

        # p-value calculation
        ods <- computePvalues(ods, alternative = "two.sided", method="BY", BPPARAM=MulticoreParam(workers = n_parallel_workers))

        # Z-score calculation
        ods <- computeZscores(ods)

        # Get results
        res <- results(ods, padjCutoff = pval_cutoff, zScoreCutoff = z_cutoff)
        res_all <- results(ods, all = TRUE)
        write_csv(res, file = "results/results.significant.csv")
        write_csv(res_all, file = "results/results.all.csv")

        # Get number of aberrant genes per sample and per gene
        ab_genes_per_samples <- as.data.frame(aberrant(ods, by = "sample", padjCutoff = 0.05, zScoreCutoff = 0))
        colnames(ab_genes_per_samples) <- "num_aberrant_genes"
        ab_genes_per_samples\\$sampleID <- rownames(ab_genes_per_samples)
        ab_genes_per_samples <- ab_genes_per_samples[, c(2, 1)]
        write_csv(
            ab_genes_per_samples,
            file = "results/aberrant_genes_per_sample.csv"
        )

        ab_samples_per_gene <- as.data.frame(aberrant(ods, by = "gene", padjCutoff = 0.05, zScoreCutoff = 0))
        colnames(ab_samples_per_gene) <- "num_aberrant_samples"
        ab_samples_per_gene\\$geneID <- rownames(ab_samples_per_gene)
        ab_samples_per_gene <- ab_samples_per_gene[, c(2, 1)]
        write_csv(
            ab_samples_per_gene,
            file = "results/aberrant_samples_per_gene.csv"
        )

        png(file = "plots/stats/aberrant_genes_per_sample.png", width = 4000, height = 4000, res = 600)
        plotAberrantPerSample(ods, padjCutoff = 0.05)
        dev.off()

        # Volcano plots
        sig_samples <- unique(res\\$sampleID)
        for (sig_sample in sig_samples) {
            png(file = paste0("plots/volcano/volcano.", sig_sample, ".png"), width = 4000, height = 4000, res = 600)
            p <- plotVolcano(ods, sig_sample, basePlot = TRUE)
            print(p)
            dev.off()
        }

        # Gene-level plots
        sig_genes <- unique(res\\$geneID)
        for (sig_gene in sig_genes) {
            png(file = paste0("plots/sig_genes/expression_rank.", sig_gene, ".png"), width = 4000, height = 4000, res = 600)
            p <- plotExpressionRank(ods, sig_gene, basePlot = TRUE)
            print(p)
            dev.off()
            png(file = paste0("plots/sig_genes/qq.", sig_gene, ".png"), width = 4000, height = 4000, res = 600)
            p <- plotQQ(ods, sig_gene)
            print(p)
            dev.off()
            png(file = paste0("plots/sig_genes/expected_vs_observed_counts.", sig_gene, ".png"), width = 4000, height = 4000, res = 600)
            p <- plotExpectedVsObservedCounts(ods, sig_gene, basePlot = TRUE)
            print(p)
            dev.off()
        }

        # Power analysis plot
        png(file = "plots/stats/power_analysis.png", width = 4000, height = 4000, res = 600)
        plotPowerAnalysis(ods)
        dev.off()

        # Save data
        save(
            ods_orig,
            ods_filt,
            ods_subset,
            ods_opt_q,
            ods_controlled,
            ods,
            res,
            res_all,
            ab_genes_per_samples,
            ab_samples_per_gene,
            file = "outrider.results.RData"
        )
        EOF
        """
        # Copy outputs
        self.command += f"""
        tar -czvf {self.output['heatmaps.tar.gz']} -C plots/heatmaps .
        tar -czvf {self.output['volcano_plots.tar.gz']} -C plots/volcano .
        tar -czvf {self.output['gene_plots.tar.gz']} -C plots/sig_genes .
        tar -czvf {self.output['stats_plots.tar.gz']} -C plots/stats .
        cp results/results.significant.csv {self.output['results.txt']}
        cp results/results.all.csv {self.output['results.all.txt']}
        cp results/aberrant_genes_per_sample.csv {self.output['aberrant_genes_per_sample.txt']}
        cp results/aberrant_samples_per_gene.csv {self.output['aberrant_samples_per_gene.txt']}
        cp outrider.results.RData {self.output['RData']}
        """
        self.command = dedent(self.command).strip()

    def __str__(self):
        return self.command

    def __repr__(self):
        return self.__str__()


def outrider(
    b: hb.Batch,
    input_counts: list[str | Path],
    output_rdata_path: str | Path | None = None,
    cohort_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> Job | None:
    """
    Run Outrider.
    """
    # Reuse existing output if possible
    if output_rdata_path and can_reuse(output_rdata_path, overwrite):
        return None

    # Localise input files
    assert all([isinstance(f, (str, Path)) for f in input_counts])
    infiles = {basename(str(f)).replace('.count', ''): str(f) for f in input_counts}
    infiles_rg = b.read_input_group(**infiles)
    infiles_localised = [str(infiles_rg[key]) for key in infiles.keys()]
    gtf_file = get_config()['references']['star'].get('gtf')
    gtf_file = to_path(gtf_file)
    gtf_file_rg = b.read_input_group(gtf=str(gtf_file))

    # Get Outrider parameters
    pval_cutoff = get_config().get('outrider', {}).get('pval_cutoff', 0.05)
    z_cutoff = get_config().get('outrider', {}).get('z_cutoff', 0.0)

    # Create job
    job_name = f'outrider_{cohort_name}' if cohort_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='outrider')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('outrider'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(j, ncpu=nthreads, storage_gb=50)

    j.declare_resource_group(
        output={
            'results.txt': '{root}.results.txt',
            'results.all.txt': '{root}.results.all.txt',
            'aberrant_genes_per_sample.txt': '{root}.aberrant_genes_per_sample.txt',
            'aberrant_samples_per_gene.txt': '{root}.aberrant_samples_per_gene.txt',
            'RData': '{root}.RData',
            'heatmaps.tar.gz': '{root}.heatmaps.tar.gz',
            'volcano_plots.tar.gz': '{root}.volcano_plots.tar.gz',
            'gene_plots.tar.gz': '{root}.gene_plots.tar.gz',
            'stats_plots.tar.gz': '{root}.stats_plots.tar.gz',
        },
    )

    # Create counting command
    outrider = Outrider(
        input_counts=infiles_localised,
        gtf_file=str(gtf_file_rg.gtf),
        output=j.output,
        nthreads=res.get_nthreads(),
        pval_cutoff=pval_cutoff,
        z_cutoff=z_cutoff,
    )

    cmd = str(outrider)
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_rdata_path:
        # NOTE: j.output is just a placeholder
        b.write_output(j.output, str(to_path(output_rdata_path).with_suffix('')))

    return j
