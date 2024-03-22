"""
Perform aberrant splicing analysis with FRASER.
"""

from textwrap import dedent

import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config, image_path
from cpg_utils.hail_batch import command
from cpg_workflows.filetypes import (
    BamPath,
    CramPath,
)
from cpg_workflows.jobs.bam_to_cram import cram_to_bam
from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse
from cpg_workflows.workflow import (
    SequencingGroup,
)


class Fraser:
    """
    Construct a FRASER command for performing aberrant splicing analysis.
    """

    def __init__(
        self,
        fds_tar: hb.ResourceFile,
        cohort_name: str,
        output: hb.ResourceGroup,
        nthreads: int = 8,
        pval_cutoff: float = 0.05,
        z_cutoff: float | None = None,
        min_delta_psi: float = 0.0,
        delta_psi_cutoff: float = 0.3,
        min_count: int = 5,
    ) -> None:
        self.fds_tar = fds_tar
        assert isinstance(self.fds_tar, hb.ResourceFile), f'fds_tar must be a resource file, instead got {self.fds_tar}'
        self.cohort_name = cohort_name
        self.output = output
        self.nthreads = nthreads
        self.pval_cutoff = str(pval_cutoff)
        self.z_cutoff = str(z_cutoff) if z_cutoff else 'NA'
        self.min_delta_psi = str(min_delta_psi)
        self.delta_psi_cutoff = str(delta_psi_cutoff)
        self.min_count = str(min_count)

        # Build OUTRIDER command
        self.command = f"""\
        # Create output directories
        rm -rf plots results output
        mkdir -p plots/heatmaps plots/volcano plots/misc results

        # Extract FDS data
        tar -xvf {self.fds_tar}
        """
        self.command += """\
        R --vanilla <<EOF
        library(FRASER)
        library(tidyverse)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        library(org.Hs.eg.db)
        """
        self.command += f"""\
        # Set significance values and other parameters
        pval_cutoff <- {self.pval_cutoff}
        z_cutoff <- {self.z_cutoff}
        minDeltaPsi <- {self.min_delta_psi}
        deltaPsi_cutoff <- {self.delta_psi_cutoff}
        min_count <- {self.min_count}
        n_parallel_workers <- {str(self.nthreads - 1)}

        # Load FDS (pre-counted)
        fds <- loadFraserDataSet(dir = "output", name = "{self.cohort_name}")

        # Extract count data and build new FDS using the FraserDataSet command
        sample_table <- fds@colData
        raw_counts_splice_sites <- cbind(
          as.data.table(granges(rowRanges(fds, type="ss"))),
          as.data.table(rowData(fds, type="ss")),
          as.data.table(counts(fds, type="ss"))
        )
        raw_counts_junctions <- cbind(
          as.data.table(granges(rowRanges(fds, type="j"))),
          as.data.table(counts(fds, type="j")),
          as.data.table(rowData(fds, type="j"))
        )
        rm(fds)

        fds <- FraserDataSet(
          colData = sample_table,
          junctions = raw_counts_junctions,
          spliceSites = raw_counts_splice_sites,
          workingDir = "output.cohort.new"
        )
        """
        self.command += """\
        # Setup parallelisation
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        # Calculate PSI values
        fds <- calculatePSIValues(fds)

        # Filter
        fds <- filterExpressionAndVariability(fds, minDeltaPsi = minDeltaPsi, filter = FALSE)

        png(file = "plots/misc/filter_expression.png", width = 4000, height = 4000, res = 600)
        plotFilterExpression(fds, bins = 100)
        dev.off()

        fds_filtered <- fds[mcols(fds, type = "j")[, "passed"],]

        # Plot sample correlation
        for (psi_type in c("psi5", "psi3", "theta")) {
            png(file = paste0("plots/heatmaps/count_correlation_heatmap.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
            print(plotCountCorHeatmap(fds_filtered, type = psi_type, logit = TRUE, normalized = FALSE))
            dev.off()
        }

        # Calculate optimal encoding dimensions (q)
        for (psi_type in c("psi5", "psi3", "theta")) {
            fds_filtered <- optimHyperParams(fds_filtered, type = psi_type, plot = FALSE)
            png(file = paste0("plots/misc/optimal_q.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
            print(plotEncDimSearch(fds_filtered, type = psi_type))
            dev.off()
        }

        optimal_qs <- c(
          psi5 = fds_filtered@metadata\\$hyperParams_psi5\\$q,
          psi3 = fds_filtered@metadata\\$hyperParams_psi3\\$q,
          theta = fds_filtered@metadata\\$hyperParams_theta\\$q
        )

        # Fit model
        fds_filtered_fit <- FRASER(fds_filtered, q = optimal_qs, BPPARAM = bp)

        # Plot sample correlation post-correction
        for (psi_type in c("psi5", "psi3", "theta")) {
            png(file = paste0("plots/heatmaps/count_correlation_heatmap.corrected.", psi_type, ".png"), width = 4000, height = 4000, res = 600)
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
            data.frame(res),
            file = paste0(
                "results/results.significant.p_", pval_cutoff, ".z_", z_cutoff, ".dPsi_", deltaPsi_cutoff, ".min_count_", min_count, ".csv"
            )
        )
        write_csv(data.frame(res_all), file = "results/results.all.csv")

        # Plot results
        for (sample_id in fds\\$sampleID) {
            png(file = paste0("plots/volcano/", sample_id, ".psi5.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "psi5")
            dev.off()
            png(file = paste0("plots/volcano/", sample_id, ".psi3.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "psi3")
            dev.off()
            png(file = paste0("plots/volcano/", sample_id, ".theta.png"), width = 4000, height = 4000, res = 600)
            plotVolcano(fds_filtered_fit, sample_id, type = "theta")
            dev.off()
        }

        # Save
        saveFraserDataSet(fds_filtered_fit, dir = getwd(), name = "FraserDataSet")
        EOF
        """
        # Tar up outputs
        self.command += f"""\
        tar -czvf {self.output['heatmaps.tar.gz']} -C plots/heatmaps .
        tar -czvf {self.output['volcano_plots.tar.gz']} -C plots/volcano .
        tar -czvf {self.output['misc_plots.tar.gz']} -C plots/misc .
        cp results/results.significant.p_{self.pval_cutoff}.z_{self.z_cutoff}.dPsi_{self.delta_psi_cutoff}.min_count_{self.min_count}.csv {self.output['results.csv']}
        cp results/results.all.csv {self.output['results.all.csv']}
        tar -czvf {self.output['fds.tar.gz']} savedObjects/
        """
        self.command = dedent(self.command).strip()

    def __str__(self):
        return self.command

    def __repr__(self):
        return self.__str__()


def fraser(
    b: hb.Batch,
    input_bams_or_crams: list[tuple[BamPath, None] | tuple[CramPath, Path]],
    cohort_name: str,
    output_fds_path: str | Path | None = None,
    job_attrs: dict[str, str] | None = None,
    overwrite: bool = False,
    requested_nthreads: int | None = None,
) -> list[Job]:
    """
    Run FRASER.
    """
    # Reuse existing output if possible
    if output_fds_path and can_reuse(output_fds_path, overwrite):
        return []

    jobs: list[Job] = []

    # Convert CRAMs to BAMs if necessary
    input_bams_localised: dict[str, hb.ResourceFile] = {}
    for input_bam_or_cram_tuple in input_bams_or_crams:
        input_bam_or_cram = input_bam_or_cram_tuple[0]
        potential_bam_path = input_bam_or_cram_tuple[1]
        if isinstance(input_bam_or_cram, CramPath):
            sample_id = input_bam_or_cram.path.name.replace('.cram', '')
            output_bam_path: Path | None = None
            if potential_bam_path:
                output_bam_path = potential_bam_path
            j, output_bam = cram_to_bam(
                b=b,
                input_cram=input_bam_or_cram.resource_group(b),
                output_bam=output_bam_path,
                job_attrs=job_attrs,
                requested_nthreads=requested_nthreads,
            )
            if j and isinstance(j, Job):
                jobs.append(j)
            input_bams_localised[sample_id] = output_bam.bam
        elif isinstance(input_bam_or_cram, BamPath):
            sample_id = input_bam_or_cram.path.name.replace('.bam', '')
            # Localise BAM
            input_bams_localised[sample_id] = input_bam_or_cram.resource_group(b).bam
    assert all([isinstance(f, hb.ResourceFile) for f in list(input_bams_localised.values())])

    # Create FRASER job
    job_name = f'fraser_{cohort_name}' if cohort_name else 'fraser'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    j.declare_resource_group(
        output={
            'results.csv': '{root}.results.csv',
            'results.all.csv': '{root}.results.all.csv',
            'heatmaps.tar.gz': '{root}.heatmaps.tar.gz',
            'volcano_plots.tar.gz': '{root}.volcano_plots.tar.gz',
            'misc_plots.tar.gz': '{root}.misc_plots.tar.gz',
            'fds.tar.gz': '{root}.fds.tar.gz',
        },
    )

    # Perform counting in parallel jobs
    output_counts_prefix = to_path(output_fds_path).parent / 'counts'
    count_jobs, fds_tar = fraser_count(
        b=b,
        input_bams_localised=input_bams_localised,
        cohort_name=cohort_name,
        output_counts_prefix=output_counts_prefix,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )
    jobs.extend(count_jobs)

    # Get Fraser parameters
    pval_cutoff = get_config().get('fraser', {}).get('pval_cutoff', 0.05)
    z_cutoff = get_config().get('fraser', {}).get('z_cutoff', None)
    min_delta_psi = get_config().get('fraser', {}).get('min_delta_psi', 0.0)
    delta_psi_cutoff = get_config().get('fraser', {}).get('delta_psi_cutoff', 0.3)
    min_count = get_config().get('fraser', {}).get('min_count', 5)

    # Create counting command
    fraser = Fraser(
        fds_tar=fds_tar,
        cohort_name=cohort_name,
        output=j.output,
        nthreads=res.get_nthreads(),
        pval_cutoff=pval_cutoff,
        z_cutoff=z_cutoff,
        min_delta_psi=min_delta_psi,
        delta_psi_cutoff=delta_psi_cutoff,
        min_count=min_count,
    )

    cmd = str(fraser)
    j.command(command(cmd, monitor_space=True))

    jobs.append(j)

    # Write output to file
    if output_fds_path:
        # NOTE: j.output is just a placeholder
        b.write_output(
            j.output,
            str(to_path(output_fds_path).with_suffix('').with_suffix('').with_suffix('')),
        )  # Remove .fds.tar.gz suffix

    return jobs


def fraser_count(
    b: hb.Batch,
    input_bams_localised: dict[str, hb.ResourceFile],
    cohort_name: str,
    output_counts_prefix: Path,
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[list[Job], hb.ResourceFile]:
    """
    Run FRASER counting.
    """
    jobs = []
    j, fds, sample_ids = fraser_init(
        b=b,
        input_bams_localised=input_bams_localised,
        cohort_name=cohort_name,
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )
    jobs.append(j)

    split_counts_dict = {}
    for sample_id in sample_ids:
        output_counts_path = output_counts_prefix / 'output/cache/splitCounts' / f'splitCounts-{sample_id}.RDS'
        sample_j, split_counts = fraser_count_split_reads_one_sample(
            b=b,
            fds=fds,
            sample_id=sample_id,
            bam=input_bams_localised[sample_id],
            cohort_name=cohort_name,
            output_counts_path=output_counts_path,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        if sample_j:
            jobs.append(sample_j)
        split_counts_dict[sample_id] = split_counts

    j, split_counts_rg = fraser_merge_split_reads(
        b=b,
        fds=fds,
        cohort_name=cohort_name,
        split_counts_dict=split_counts_dict,
        bams=list(input_bams_localised.values()),
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )
    jobs.append(j)

    non_spliced_counts_dict = {}
    for sample_id in sample_ids:
        output_counts_path = (
            output_counts_prefix / 'output/cache/nonSplicedCounts' / cohort_name / f'nonSplicedCounts-{sample_id}.h5'
        )
        sample_j, non_spliced_counts = fraser_count_non_split_reads_one_sample(
            b=b,
            fds=fds,
            cohort_name=cohort_name,
            split_counts=split_counts_rg,
            sample_id=sample_id,
            bam=input_bams_localised[sample_id],
            output_counts_path=output_counts_path,
            job_attrs=job_attrs,
            requested_nthreads=requested_nthreads,
        )
        jobs.append(sample_j)
        non_spliced_counts_dict[sample_id] = non_spliced_counts

    j, fds_tar = fraser_merge_non_split_reads(
        b=b,
        fds=fds,
        cohort_name=cohort_name,
        non_spliced_counts_dict=non_spliced_counts_dict,
        split_counts=split_counts_rg,
        bams=list(input_bams_localised.values()),
        job_attrs=job_attrs,
        requested_nthreads=requested_nthreads,
    )
    jobs.append(j)

    return jobs, fds_tar


def fraser_init(
    b: hb.Batch,
    input_bams_localised: dict[str, hb.ResourceFile],
    cohort_name: str,
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceFile, list[str]]:
    """
    Run FRASER initialisation.
    """
    # Create FRASER job
    job_name = 'fraser_init'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    bam_files_r_str = ''
    sample_ids_r_str = ''
    sample_ids = []
    for sample_id, bam_file in input_bams_localised.items():
        bam_files_r_str += f'"{bam_file}", '
        sample_ids_r_str += f'"{sample_id}", '
        sample_ids.append(sample_id)
    bam_files_r_str = bam_files_r_str[:-2]  # Remove trailing comma and space
    sample_ids_r_str = sample_ids_r_str[:-2]  # Remove trailing comma and space

    cmd = dedent(
        f"""\
        R --vanilla <<EOF
        library(FRASER)

        bam_files <- c({bam_files_r_str})
        sample_ids <- c({sample_ids_r_str})
        sample_table <- DataFrame(
            sampleID = sample_ids,
            bamFile = bam_files,
            group = seq_len(length(bam_files)),
            pairedEnd = TRUE
        )

        fds <- FraserDataSet(
            colData = sample_table,
            workingDir = "output",
            name = "{cohort_name}"
        )

        n_parallel_workers <- {str(res.get_nthreads() - 1)}
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        fds <- saveFraserDataSet(fds)
        EOF

        # Move output to resource file
        mv output/savedObjects/{cohort_name}/fds-object.RDS {j.fds}
        """,
    )

    j.command(command(cmd, monitor_space=True))

    return j, j.fds, sample_ids


def fraser_count_split_reads_one_sample(
    b: hb.Batch,
    fds: hb.ResourceFile,
    sample_id: str,
    bam: hb.ResourceFile,
    cohort_name: str,
    output_counts_path: Path,
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job | None, hb.ResourceFile]:
    """
    Run FRASER split-read counting for a single sample.
    """
    # Reuse existing output if possible
    if can_reuse(output_counts_path):
        return None, b.read_input(str(output_counts_path))

    # Create FRASER job
    job_name = 'fraser_count_split'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    cmd = dedent(
        f"""\
        # Symlink FDS resource file to proper location
        mkdir -p output/savedObjects/{cohort_name}
        ln -s {fds} output/savedObjects/{cohort_name}/fds-object.RDS
        # ls BAM file to ensure it is localised
        ls {bam}

        R --vanilla <<EOF
        library(FRASER)

        fds <- loadFraserDataSet(dir = "output", name = "{cohort_name}")

        n_parallel_workers <- {str(res.get_nthreads() - 1)}
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        options("FRASER.maxSamplesNoHDF5"=0)
        options("FRASER.maxJunctionsNoHDF5"=-1)

        sample_id <- "{sample_id}"
        sample_count <- countSplitReads(
          sampleID = sample_id,
          fds = fds,
          NcpuPerSample = n_parallel_workers,
        )
        EOF

        # Move output to resource file
        mv output/cache/splitCounts/splitCounts-{sample_id}.RDS {j.split_counts}
        """,
    )

    j.command(command(cmd, monitor_space=True))

    # Write output to file
    b.write_output(j.split_counts, str(output_counts_path))

    return j, j.split_counts


def fraser_merge_split_reads(
    b: hb.Batch,
    fds: hb.ResourceFile,
    cohort_name: str,
    split_counts_dict: dict[str, hb.ResourceFile],
    bams: list[hb.ResourceFile],
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceGroup]:
    """
    Merge split-read counts.
    """
    # Create FRASER job
    job_name = 'fraser_merge_split'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    # Create command to symlink split counts
    link_counts_cmd = 'mkdir -p output/cache/splitCounts\n'
    link_counts_cmd += '\n'.join(
        [
            f'ln -s {split_counts} output/cache/splitCounts/splitCounts-{sample_id}.RDS'
            for sample_id, split_counts in split_counts_dict.items()
        ],
    )

    # Create resource group for outputs
    split_counts_rg = {
        'raw_counts_j_h5': f'output/savedObjects/{cohort_name}/rawCountsJ.h5',
        'split_counts_assays': f'output/savedObjects/{cohort_name}/splitCounts/assays.h5',
        'split_counts_se': f'output/savedObjects/{cohort_name}/splitCounts/se.rds',
        'g_ranges_split_counts': 'rds/g_ranges_split_counts.RDS',
        'g_ranges_non_split_counts': 'rds/g_ranges_non_split_counts.RDS',
        'splice_site_coords': 'rds/splice_site_coords.RDS',
    }
    split_counts_rg_flat = {key: value.split('/')[-1] for key, value in split_counts_rg.items()}
    j.declare_resource_group(
        split_counts=split_counts_rg_flat,
    )

    # Create move command for outputs
    move_cmd = '\n'.join([f'mv {file} {j.split_counts[key]}' for key, file in split_counts_rg.items()])

    cmd = dedent(
        f"""\
        # Symlink FDS resource file to proper location
        mkdir -p output/savedObjects/{cohort_name}
        ln -s {fds} output/savedObjects/{cohort_name}/fds-object.RDS
        # Symlink split counts
        {link_counts_cmd}
        # ls BAM files to ensure they are localised
        ls {' '.join(bams)}
        # Make RDS directory
        mkdir -p rds

        R --vanilla <<EOF
        library(FRASER)
        fds <- loadFraserDataSet(dir = "output", name = "{cohort_name}")

        n_parallel_workers <- {str(res.get_nthreads() - 1)}
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        options("FRASER.maxSamplesNoHDF5"=0)
        options("FRASER.maxJunctionsNoHDF5"=-1)

        minExpressionInOneSample <- 20

        # Read counts from cache
        split_counts <- getSplitReadCountsForAllSamples(fds = fds, recount = FALSE)

        split_count_ranges <- rowRanges(split_counts)
        split_count_ranges <- FRASER:::annotateSpliceSite(split_count_ranges)
        saveRDS(split_count_ranges, "rds/g_ranges_split_counts.RDS")

        max_count <- rowMaxs(assay(split_counts, "rawCountsJ"))
        passed <- max_count >= minExpressionInOneSample
        split_count_ranges <- split_count_ranges[passed, ]
        saveRDS(split_count_ranges, "rds/g_ranges_non_split_counts.RDS")

        splice_site_coords <- FRASER:::extractSpliceSiteCoordinates(split_count_ranges, fds)
        saveRDS(splice_site_coords, "rds/splice_site_coords.RDS")
        EOF

        # Move outputs to resource group
        {move_cmd}
        """,
    )

    j.command(command(cmd, monitor_space=True))

    return j, j.split_counts


def fraser_count_non_split_reads_one_sample(
    b: hb.Batch,
    fds: hb.ResourceFile,
    cohort_name: str,
    split_counts: hb.ResourceGroup,
    sample_id: str,
    bam: hb.ResourceFile,
    output_counts_path: Path,
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceFile]:
    """
    Run FRASER non-split-read counting for a single sample.
    """
    # NOTE: Can't reuse output, need to count non-split reads for each sample

    # Create FRASER job
    job_name = 'fraser_count_non_split'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    cmd = dedent(
        f"""\
        # Symlink FDS resource file to proper location
        mkdir -p output/savedObjects/{cohort_name}
        ln -s {fds} output/savedObjects/{cohort_name}/fds-object.RDS
        # Symlink splice site coords RDS file
        mkdir -p rds
        ln -s {split_counts.splice_site_coords} rds/splice_site_coords.RDS
        # ls BAM file to ensure it is localised
        ls {bam}

        R --vanilla <<EOF
        library(FRASER)
        fds <- loadFraserDataSet(dir = "output", name = "{cohort_name}")

        n_parallel_workers <- {str(res.get_nthreads() - 1)}
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        options("FRASER.maxSamplesNoHDF5"=0)
        options("FRASER.maxJunctionsNoHDF5"=-1)

        splice_site_coords <- readRDS("rds/splice_site_coords.RDS")

        sample_id <- "{sample_id}"
        sample_count <- countNonSplicedReads(
          sampleID = sample_id,
          splitCountRanges = NULL,
          fds = fds,
          minAnchor = 5,
          NcpuPerSample = n_parallel_workers,
          spliceSiteCoords = splice_site_coords
        )
        EOF

        # Move output to resource file
        mv output/cache/nonSplicedCounts/{cohort_name}/nonSplicedCounts-{sample_id}.h5 {j.non_spliced_counts}
        """,
    )

    j.command(command(cmd, monitor_space=True))

    # Write output to file
    b.write_output(j.non_spliced_counts, str(output_counts_path))

    return j, j.non_spliced_counts


def fraser_merge_non_split_reads(
    b: hb.Batch,
    fds: hb.ResourceFile,
    cohort_name: str,
    non_spliced_counts_dict: dict[str, hb.ResourceFile],
    split_counts: hb.ResourceGroup,
    bams: list[hb.ResourceFile],
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> tuple[Job, hb.ResourceFile]:
    """
    Merge non-split-read counts.
    """
    # Create FRASER job
    job_name = 'fraser_merge_non_split'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='fraser')
    j = b.new_job(job_name, _job_attrs)
    j.image(image_path('fraser'))

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    # Create command to symlink non-spliced counts
    link_counts_cmd = f'mkdir -p output/cache/nonSplicedCounts/{cohort_name}\n'
    link_counts_cmd += '\n'.join(
        [
            f'ln -s {non_spliced_counts} output/cache/nonSplicedCounts/{cohort_name}/nonSplicedCounts-{sample_id}.h5'
            for sample_id, non_spliced_counts in non_spliced_counts_dict.items()
        ],
    )
    link_counts_cmd += '\n'
    # Add command to symlink RDS files
    link_counts_cmd += 'mkdir -p rds\n'
    link_counts_cmd += f'ln -s {split_counts.g_ranges_split_counts} rds/g_ranges_split_counts.RDS\n'
    link_counts_cmd += f'ln -s {split_counts.g_ranges_non_split_counts} rds/g_ranges_non_split_counts.RDS\n'
    link_counts_cmd += f'ln -s {split_counts.splice_site_coords} rds/splice_site_coords.RDS\n'
    # Add command to symlink split counts assays
    link_counts_cmd += f'mkdir -p output/savedObjects/{cohort_name}/splitCounts\n'
    link_counts_cmd += f'ln -s {split_counts.raw_counts_j_h5} output/savedObjects/{cohort_name}/rawCountsJ.h5\n'
    link_counts_cmd += (
        f'ln -s {split_counts.split_counts_assays} output/savedObjects/{cohort_name}/splitCounts/assays.h5\n'
    )
    link_counts_cmd += f'ln -s {split_counts.split_counts_se} output/savedObjects/{cohort_name}/splitCounts/se.rds\n'

    cmd = dedent(
        f"""\
        # Symlink FDS resource file to proper location
        mkdir -p output/savedObjects/{cohort_name}
        ln -s {fds} output/savedObjects/{cohort_name}/fds-object.RDS
        # Symlink non-spliced counts, RDS files, and split counts
        {link_counts_cmd}
        # ls BAM files to ensure they are localised
        ls {' '.join(bams)}

        R --vanilla <<EOF
        library(FRASER)

        fds <- loadFraserDataSet(dir = "output", name = "{cohort_name}")

        n_parallel_workers <- {str(res.get_nthreads() - 1)}
        register(MulticoreParam(workers = n_parallel_workers))
        bp <- MulticoreParam(workers = n_parallel_workers)

        options("FRASER.maxSamplesNoHDF5"=0)
        options("FRASER.maxJunctionsNoHDF5"=-1)

        minAnchor <- 5

        split_count_ranges <- readRDS("rds/g_ranges_non_split_counts.RDS")
        non_split_counts <- getNonSplitReadCountsForAllSamples(
          fds = fds,
          splitCountRanges = split_count_ranges,
          minAnchor = minAnchor,
          recount = FALSE,
          longRead = FALSE
        )

        split_count_ranges <- readRDS("rds/g_ranges_split_counts.RDS")
        splice_site_coords <- readRDS("rds/splice_site_coords.RDS")

        split_counts_h5 <- HDF5Array::HDF5Array("output/savedObjects/{cohort_name}/rawCountsJ.h5", "rawCountsJ")
        split_counts_se <- SummarizedExperiment(
          colData = colData(fds),
          rowRanges = split_count_ranges,
          assays = list(rawCountsJ = split_counts_h5)
        )

        non_split_counts_h5 <- HDF5Array::HDF5Array("output/savedObjects/{cohort_name}/rawCountsSS.h5", "rawCountsSS")
        non_split_counts_se <- SummarizedExperiment(
          colData = colData(fds),
          rowRanges = splice_site_coords,
          assays = list(rawCountsSS = non_split_counts_h5)
        )

        fds <- addCountsToFraserDataSet(
          fds = fds,
          splitCounts = split_counts_se,
          nonSplitCounts = non_split_counts_se
        )

        fds <- saveFraserDataSet(fds)
        EOF

        # tar saved objects directory
        tar -chf {j.fds_tar} output/savedObjects/{cohort_name}/
        """,
    )

    j.command(command(cmd, monitor_space=True))

    return j, j.fds_tar
