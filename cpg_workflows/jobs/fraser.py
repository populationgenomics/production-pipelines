"""
Perform aberrant splicing analysis with FRASER.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command, image_path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
from cpg_workflows.jobs.bam_to_cram import cram_to_bam
from textwrap import dedent
from os.path import basename


class Fraser:
    """
    Construct a FRASER command for performing aberrant splicing analysis.
    """

    def __init__(
            self,
            input_bams: list[str | Path],
            output_tar_gz_path: str | Path | None = None,
            nthreads: int = 8,
        ) -> None:
        self.input_bams = input_bams
        assert isinstance(self.input_bams, list), f'input_bams must be a list, instead got {self.input_bams}'
        self.input_bams_r_str = ', '.join([f'"{str(f)}"' for f in self.input_bams])
        self.output_tar_gz_path = str(output_tar_gz_path)
        self.nthreads = nthreads

        # Build OUTRIDER command
        self.command = """
        R --vanilla <<EOF
        library(FRASER)
        library(tidyverse)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        library(org.Hs.eg.db)

        # Create directories
        dir.create("plots", showWarnings = FALSE)
        dir.create("plots/heatmaps", showWarnings = FALSE)
        dir.create("plots/by_gene", showWarnings = FALSE)
        dir.create("plots/by_sample", showWarnings = FALSE)
        dir.create("results", showWarnings = FALSE)
        dir.create("output", showWarnings = FALSE)
        """
        self.command += f"""
        # Set significance values
        pval_cutoff <- 0.05
        z_cutoff <- NA
        minDeltaPsi <- 0.0
        deltaPsi_cutoff <- 0.3
        min_count <- 5
        n_parallel_workers <- {str(self.nthreads - 1)}

        input_bam_files <- c({self.input_bams_r_str})
        """
        self.command += """
        # Create sample table
        sampleTable <- DataFrame(
        sampleID = gsub("\\\\\\\\.bam$", "", basename(input_bam_files), perl = TRUE),
        bamFile = input_bam_files,
        group = seq_len(length(input_bam_files)),
        pairedEnd = TRUE
        )

        # Create FRASER settings object
        settings <- FraserDataSet(
        colData = sampleTable,
        workingDir = "output"
        )

        # Setup parallelisation
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        # Count reads
        fds <- countRNAData(settings)

        # Calculate PSI values
        fds <- calculatePSIValues(fds)

        # Filter
        fds <- filterExpressionAndVariability(fds, minDeltaPsi = minDeltaPsi, filter = FALSE)

        png(file = "plots/filter_expression.png", width = 4000, height = 4000, res = 600)
        plotFilterExpression(fds, bins = 100)
        dev.off()

        fds_filtered <- fds[mcols(fds, type = "j")[, "passed"],]

        # Plot sample correlation
        for (psi_type in c("psi5", "psi3", "theta")) {
            png(file = paste0("plots/count_correlation_heatmap.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
            print(plotCountCorHeatmap(fds_filtered, type = psi_type, logit = TRUE, normalized = FALSE))
            dev.off()
        }

        # Calculate optimal encoding dimensions (q)
        optimal_qs <- c(
            psi5 = optimHyperParams(fds_filtered, type = "psi5", plot = FALSE),
            psi3 = optimHyperParams(fds_filtered, type = "psi3", plot = FALSE),
            theta = optimHyperParams(fds_filtered, type = "theta", plot = FALSE)
        )
        for (psi_type in c("psi5", "psi3", "theta")) {
            png(file = paste0("plots/optimal_q.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
            print(plotEncDimSearch(fds_filtered, type = psi_type))
            dev.off()
        }

        # Fit model
        fds_filtered_fit <- FRASER(fds_filtered, q = optimal_qs, BPPARAM = bp)

        # Plot sample correlation post-correction
        for (psi_type in c("psi5", "psi3", "theta")) {
            png(file = paste0("plots/count_correlation_heatmap.corrected.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
            print(plotCountCorHeatmap(fds_filtered_fit, type = psi_type, logit = TRUE, normalized = TRUE))
            dev.off()
        }

        # Annotate gene symbols
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        orgDb <- org.Hs.eg.db
        fds_filtered_fit <- annotateRangesWithTxDb(fds_filtered_fit, txdb = txdb, orgDb = orgDb)

        # Get results
        res <- results(fds_filtered_fit, padjCutoff = pval_cutoff, deltaPsiCutoff = deltaPsi_cutoff, zScoreCutoff = z_cutoff, minCount = min_count)
        res_all <- results(fds_filtered_fit, padjCutoff = 1, deltaPsiCutoff = 0, minCount = 0)
        write_csv(
            res,
            file = paste0(
                "results/results.significant.p_", pval_cutoff, ".z_", z_cutoff, ".dPsi_", deltaPsi_cutoff, ".min_count_", min_count, ".csv"
            )
        )
        write_csv(res_all, file = "results/results.all.csv")

        # Plot results
        for (sample_id in sampleTable\\$sampleID) {
            png(file = paste0("plots/per_sample/", sample_id, ".psi5.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "psi5")
            dev.off()
            png(file = paste0("plots/per_sample/", sample_id, ".psi3.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "psi3")
            dev.off()
            png(file = paste0("plots/per_sample/", sample_id, ".theta.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "theta")
            dev.off()
        }

        # Save
        saveFraserDataSet(fds_filtered_fit, dir = getwd(), name = "FraserDataSet")
        EOF
        """
        # Tar up outputs
        self.command += f"""
        tar -czvf {self.output_tar_gz_path} plots results output savedObjects
        """
        self.command = dedent(self.command).strip()

    def __str__(self):
        return self.command
    
    def __repr__(self):
        return self.__str__()


def fraser(
    b: hb.Batch,
    input_bams_or_crams: list[BamPath | CramPath],
    output_path: str | Path | None = None,
    cohort_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> list[Job]:
    """
    Run FRASER.
    """
    # Reuse existing output if possible
    if output_path and can_reuse(output_path, overwrite):
        return []

    jobs: list[Job] = []

    # Convert CRAMs to BAMs if necessary
    input_bams_localised: list[hb.ResourceGroup] = []
    for input_bam_or_cram in input_bams_or_crams:
        if isinstance(input_bam_or_cram, CramPath):
            j, output_bam = cram_to_bam(
                b=b,
                input_cram=input_bam_or_cram,
                job_attrs=job_attrs,
                requested_nthreads=requested_nthreads,
            )
            if j and isinstance(j, Job):
                jobs.append(j)
            input_bam_or_cram = output_bam
        elif isinstance(input_bam_or_cram, BamPath):
            # Localise BAM
            input_bams_localised.append(input_bam_or_cram.resource_group(b))
    assert all([isinstance(f, hb.ResourceGroup) for f in input_bams_localised])

    # Create FRASER job
    job_name = f'fraser_{cohort_name}' if cohort_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    # j.image(image_path('fraser'))
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/fraser:1.12.1')

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    j.declare_resource_group(
        output={
            'tar_gz': '{root}.tar.gz'
        }
    )

    # Create counting command
    fraser = Fraser(
        input_bams=input_bams_localised,
        output_tar_gz_path=j.output.tar_gz,
        nthreads=res.get_nthreads(),
    )

    cmd = str(fraser)
    j.command(command(cmd, monitor_space=True))

    jobs.append(j)

    # Write output to file
    if output_path:
        # NOTE: j.output is just a placeholder
        b.write_output(j.output.tar_gz, str(output_path))
    
    return jobs
