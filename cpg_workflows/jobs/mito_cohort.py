"""
Create Hail Batch jobs to call mitochondrial SNVs
"""
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils.config import get_config
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse

from cpg_workflows.mito_pipeline_scripts import (
    annotate_coverage as annotate_coverage_script,
)


def annotate_coverage(
    b,
    base_level_coverage_by_sgid: dict[str, Path],
    coverage_ht: Path,
    checkpoint_prefix: Path,
    # haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Implementation of Broad script here: https://github.com/broadinstitute/gnomad-mitochondria/blob/main/gnomad_mitochondria/pipeline/annotate_coverage.py

    Combines the per base coverage files from each sequence group into a mt(MatrixTable), ht(HailTable), and tsv file, and will also calculate the following aggregate statistics per base:
        Mean coverage
        Median coverage
        Fraction of samples with > 100x coverage
        Fraction of samples with > 10000x coverage

        Required inputs:
            base_level_coverage_by_sid: dict of path to tab-delimited coverage file for each sequence group.
            coverage_ht: path to write output ht

    """
    job_attrs = job_attrs or {}
    j = b.new_python_job('mito_annotate_coverage', job_attrs)
    j.image(get_config()['workflow']['driver_image'])

    res = STANDARD.request_resources(ncpu=8)
    res.set_to_job(j)

    # generate tsv string to use as input file
    # required format: "tsv of participant_id, base_level_coverage_metrics, sample"
    tsv_string = "#participant_id\\tpath\\tsample_id"
    inputs = []
    for sid, path in base_level_coverage_by_sgid.items():
        tsv_string += f'{sid}\\t{path}\\t{sid}\\n'
        inputs.append((sid, path, sid))

    # script needs to write to a .ht path
    # j.declare_resource_group(
    #     outfile={'ht': '{root}.ht'}
    # )

    ##################
    # Sticking all nested functions here...
    ##################
    import logging
    import math
    import re
    import sys
    import os

    import hail as hl

    from hail.utils.java import info

    logging.basicConfig(
        format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("Annotate coverage")
    logger.setLevel(logging.INFO)

    def multi_way_union_mts(
        mts: list, checkpoint_prefix: str, chunk_size: int
    ) -> hl.MatrixTable:
        """
        Hierarchically join together MatrixTables in the provided list.

        :param mts: List of MatrixTables to join together
        :param checkpoint_prefix: Path to temporary directory for intermediate results
        :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
        :return: Joined MatrixTable
        """
        # Convert the MatrixTables to tables where entries are an array of structs
        staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
        stage = 0
        while len(staging) > 1:
            logger.info(f"len(staging): {len(staging)}")
            # Calculate the number of jobs to run based on the chunk size
            n_jobs = int(math.ceil(len(staging) / chunk_size))
            info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
            next_stage = []

            for i in range(n_jobs):
                # Grab just the tables for the given job
                to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
                info(
                    f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
                )

                # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
                merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
                # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
                merged = merged.annotate(
                    __entries=hl.flatten(
                        hl.range(hl.len(merged.__entries)).map(
                            # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                            lambda i: hl.coalesce(
                                merged.__entries[i].__entries,
                                hl.range(hl.len(merged.__cols[i].__cols)).map(
                                    lambda j: hl.missing(
                                        merged.__entries.__entries.dtype.element_type.element_type
                                    )
                                ),
                            )
                        )
                    )
                )

                logger.info(f"merged.count(): {merged.count()}")
                # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
                merged = merged.annotate_globals(
                    __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
                )
                merged.show(10)
                next_stage.append(
                    merged.checkpoint(
                        os.path.join(checkpoint_prefix, f"stage_{stage}_job_{i}.ht"),
                        overwrite=True,
                    )
                )
            info(f"Completed stage {stage}")
            stage += 1
            staging.clear()
            staging.extend(next_stage)

        # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
        return (
            staging[0]
            ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
            .unfilter_entries()
        )

    def annotate_coverage_worker(
        inputs: list,
        output_ht: str,
        checkpoint_prefix: Path,
        chunk_size: int = 100,
        overwrite=True,
    ):
        from cpg_utils.hail_batch import init_batch

        init_batch()

        if overwrite is False and hl.hadoop_exists(output_ht):
            logger.warning(
                "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
                output_ht,
            )
        # Ensure that user supplied ht extension for output_ht
        if not output_ht.endswith(".ht"):
            sys.exit("Path supplied as output_ht must end with .ht extension")

        mt_list = []
        logger.info(
            "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
        )
        for line in inputs:
            participant_id, base_level_coverage_metrics, sample = line[0:3]
            logger.info(f"starting import of {base_level_coverage_metrics}")
            mt = hl.import_matrix_table(
                str(base_level_coverage_metrics),
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
            ).drop("target")
            mt = mt.rename({"x": "coverage"})
            mt = mt.key_cols_by(s=sample)
            mt_list.append(mt)
            logger.info(f"Finished import of {base_level_coverage_metrics}")

        logger.info("Joining individual coverage mts...")

        cov_mt = multi_way_union_mts(mt_list, str(checkpoint_prefix), chunk_size)
        n_samples = cov_mt.count_cols()

        logger.info("Adding coverage annotations...")
        # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
        cov_mt = cov_mt.annotate_rows(
            locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
            mean=hl.float(hl.agg.mean(cov_mt.coverage)),
            median=hl.median(hl.agg.collect(cov_mt.coverage)),
            over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
            over_1000=hl.float(
                (hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)
            ),
        )
        cov_mt.show()

        cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

        output_mt = re.sub(r"\.ht$", ".mt", output_ht)
        output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
        output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

        print(f"Writing sample level coverage to {output_samples}")
        sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
        sample_mt.coverage.export(str(output_samples))

        print("Writing coverage mt to {output_mt}")
        cov_mt.write(output_mt, overwrite=overwrite)

        cov_ht = cov_mt.rows()
        print("Writing coverage ht to {output_ht}")
        cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)

        print("Writing coverage tsv to {output_tsv}")
        cov_ht.export(str(output_tsv))

    # Call the worker job
    j.call(annotate_coverage_worker, inputs, str(coverage_ht), checkpoint_prefix)

    return j


def combine_vcfs(
    b,
    vcf_path_by_sgid: dict[str, Path],
    coverage_mt_path: Path,
    artifact_prone_sites_path: Path,
    combined_vcf_mt_path: Path,
    checkpoint_prefix: Path,
    chunk_size: int = 100,
    minimum_homref_coverage: int = 100,
    job_attrs: dict | None = None,
) -> Job:
    """
    Implementation of Broad script here: https://github.com/broadinstitute/gnomad-mitochondria/blob/main/gnomad_mitochondria/pipeline/combine_vcfs.py

    Takes individual sample vcfs and combines them into one vcf/mt. This script:
        - Combines individual sample vcfs into one vcf
        - Removes the "possible_numt", "mt_many_low_hets", and "FAIL" FT filters because
            these filters were found to have low performance but still may be present if
            running earlier versions of the Mutect2 mitochondria pipeline
        - For sites without a call in a particular sample, sets the call to missing if the
            sample has coverage <= minimum_homref_coverage (default: 100x) at that
            position, otherwise sets the call to homoplasmic reference
        - Applies the "artifact_prone_site" filter to any SNP or deletion that spans a
            known problematic site supplied in the "artifact_prone_sites_path" bed file

    """
    job_attrs = job_attrs or {}
    j = b.new_python_job('mito_combine_vcfs', job_attrs)
    j.image(get_config()['workflow']['driver_image'])

    res = STANDARD.request_resources(ncpu=8)
    res.set_to_job(j)

    ##################
    # python jobs seem to need called functions to be defined in this scope
    # Not sure how to get around them so jamming them all here...
    ##################
    import logging
    import math
    import os
    import re

    import hail as hl
    from hail.utils.java import info
    from typing import Dict

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(asctime)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("combine_mitochondria_vcfs_into_mt")
    logger.setLevel(logging.INFO)

    # logger.info("Setting hail flag to avoid array index out of bounds error...")
    # Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
    # TODO: reassess if this flag is still needed for future versions of Hail #CS trying without this flag...
    # hl._set_flags(no_whole_stage_codegen="1")

    def multi_way_union_mts(
        mts: list, checkpoint_prefix: str, chunk_size: int
    ) -> hl.MatrixTable:
        """
        Hierarchically join together MatrixTables in the provided list.

        :param mts: List of MatrixTables to join together
        :param checkpoint_prefix: Path to temporary directory for intermediate results
        :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
        :return: Joined MatrixTable
        """
        # Convert the MatrixTables to tables where entries are an array of structs
        staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
        stage = 0
        while len(staging) > 1:
            # Calculate the number of jobs to run based on the chunk size
            n_jobs = int(math.ceil(len(staging) / chunk_size))
            info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
            next_stage = []

            for i in range(n_jobs):
                # Grab just the tables for the given job
                to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
                info(
                    f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
                )

                # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
                merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
                # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
                merged = merged.annotate(
                    __entries=hl.flatten(
                        hl.range(hl.len(merged.__entries)).map(
                            # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                            lambda i: hl.coalesce(
                                merged.__entries[i].__entries,
                                hl.range(hl.len(merged.__cols[i].__cols)).map(
                                    lambda j: hl.null(
                                        merged.__entries.__entries.dtype.element_type.element_type
                                    )
                                ),
                            )
                        )
                    )
                )

                # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
                merged = merged.annotate_globals(
                    __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
                )

                next_stage.append(
                    merged.checkpoint(
                        os.path.join(checkpoint_prefix, f"stage_{stage}_job_{i}.ht"),
                        overwrite=True,
                    )
                )
            info(f"Completed stage {stage}")
            stage += 1
            staging.clear()
            staging.extend(next_stage)

        # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
        return (
            staging[0]
            ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
            .unfilter_entries()
        )

    def join_mitochondria_vcfs_into_mt(
        vcf_paths: Dict[str, str], checkpoint_prefix: str, chunk_size: int = 100
    ) -> hl.MatrixTable:
        """
        Reformat and join individual mitochondrial VCFs into one MatrixTable.

        :param vcf_paths: Dictionary of samples to combine (sample as key, path to VCF as value)
        :param checkpoint_prefix: Path to temporary directory for intermediate results
        :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
        :return: Joined MatrixTable of samples given in vcf_paths dictionary
        """
        mt_list = []
        for sample, vcf_path in vcf_paths.items():
            try:
                mt = hl.import_vcf(vcf_path, reference_genome="GRCh38")
            except Exception as e:
                raise ValueError(
                    f"vcf path {vcf_path} does not exist for sample {sample}"
                ) from e

            # Because the vcfs are split, there is only one AF value, although misinterpreted as an array because Number=A in VCF header
            # Second value of MMQ is the value of the mapping quality for the alternate allele
            # Add FT annotation for sample genotype filters (pull these from filters annotations of the single-sample VCFs)
            mt = mt.select_entries("DP", HL=mt.AF[0])
            mt = mt.annotate_entries(
                MQ=hl.float(mt.info["MMQ"][1]),
                TLOD=mt.info["TLOD"][0],
                FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters),
            )
            # Use GRCh37 reference as most external resources added in downstream scripts use GRCh37 contig names
            # (although note that the actual sequences of the mitochondria in both GRCh37 and GRCh38 are the same)
            mt = mt.key_rows_by(
                locus=hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
                alleles=mt.alleles,
            )
            mt = mt.key_cols_by(s=sample)
            mt = mt.select_rows()
            mt_list.append(mt)

        combined_mt = multi_way_union_mts(mt_list, checkpoint_prefix, chunk_size)

        return combined_mt

    def remove_genotype_filters(
        mt: hl.MatrixTable,
        filters_to_remove: set = {
            "possible_numt",
            "mt_many_low_hets",
            "FAIL",
            "blacklisted_site",
        },
    ) -> hl.MatrixTable:
        """
        Remove unneeded sample-level genotype filters (in FT field of the VCF) specified by the filters_to_remove parameter.

        By default, remove the 'possible_numt', 'mt_many_low_hets', and 'FAIL' filters because these filters were found to have low performance.
        Also remove the 'blacklisted_site' filter because this filter did not always behave as expected in early GATK versions. This filter can be reimplemented with the apply_mito_artifact_filter function.

        :param mt: MatrixTable containing genotype filters in the FT field of the VCF that should be removed
        :param filters_to_remove: List of genptype filters (in FT field of VCF) that should be removed from the entries
        :return: MatrixTable with specific genotype filters (in FT field of VCF) removed
        """
        mt = mt.annotate_entries(FT=mt.FT.difference(filters_to_remove))

        # If no filters exist after removing those specified above, set the FT field to PASS
        mt = mt.annotate_entries(FT=hl.if_else(hl.len(mt.FT) == 0, {"PASS"}, mt.FT))

        return mt

    def determine_hom_refs(
        mt: hl.MatrixTable, coverage_mt_path: str, minimum_homref_coverage: int = 100
    ) -> hl.MatrixTable:
        """
        Use coverage to distinguish between homref and missing sites.

        :param mt: MatrixTable from initial multi-sample merging, without homref sites determined
        :param coverage_mt_path: MatrixTable of sample level coverage at each position (per-sample and per-base; can be generated by running annotate_coverage.py)
        :param minimum_homref_coverage: Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing
        :return: MatrixTable with missing genotypes converted to homref depending on coverage
        """
        # Convert coverage to build GRCh37 to match contig names
        # Note: the mitochondrial reference genome is the same for GRCh38 and GRCh37
        coverages = hl.read_matrix_table(coverage_mt_path)
        coverages = coverages.key_rows_by(
            locus=hl.locus("MT", coverages.locus.position, reference_genome="GRCh37")
        )

        mt = mt.annotate_entries(
            DP=hl.if_else(
                hl.is_missing(mt.HL), coverages[mt.locus, mt.s].coverage, mt.DP
            )
        )

        hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)

        mt = mt.annotate_entries(
            HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
            FT=hl.if_else(hom_ref_expr, {"PASS"}, mt.FT),
            DP=hl.if_else(
                hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
                hl.null(hl.tint32),
                mt.DP,
            ),
        )

        return mt

    def apply_mito_artifact_filter(
        mt: hl.MatrixTable,
        artifact_prone_sites_path: str,
    ) -> hl.MatrixTable:
        """
        Add in artifact_prone_site filter.

        :param mt: MatrixTable to be annotated with artifact_prone_sites filter
        :param artifact_prone_sites_path: Path to BED file of artifact_prone_sites to flag in the filters column
        :return: MatrixTable with artifact_prone_sites filter
        """
        # Apply "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site
        bed = hl.import_bed(artifact_prone_sites_path)
        bed = bed.annotate(target="artifact")

        # Create a region annotation containing the interval that the variant overlaps (for SNP will be one position, but will be longer for deletions based on the length of the deletion)
        mt = mt.annotate_rows(
            region=hl.interval(
                hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
                hl.locus(
                    "MT",
                    mt.locus.position + hl.len(mt.alleles[0]) - 1,
                    reference_genome="GRCh37",
                ),
                includes_end=True,
            )
        )

        # Annotate if the start of the variant overlaps an interval in the bed file
        mt = mt.annotate_rows(
            start_overlaps=bed.index(mt.region.start, all_matches=True)
        )

        # Annotate if the end of the variant overlaps an interval in the bed file
        mt = mt.annotate_rows(end_overlaps=bed.index(mt.region.end, all_matches=True))

        # Create struct containing locus and allele (need to the check if any position of the allele overlaps an artifact-prone site, not just the locus)
        mt_temp = mt.annotate_rows(
            variant=hl.struct(locus=mt.locus, alleles=mt.alleles)
        )
        mt_temp = mt_temp.key_rows_by(mt_temp.region)

        # Need to account for cases where the start and end of the variant interval don't fall within a bed interval, but start before and after the interval (the bed interval falls completely within the variant interval)
        bed_temp = bed.annotate(
            contained_mt_alleles=mt_temp.index_rows(
                bed.interval.start, all_matches=True
            ).variant
        )

        # Explode so that each allele is on its own row and create locus and allele annotations
        bed_temp = bed_temp.explode(bed_temp.contained_mt_alleles).rename(
            {"contained_mt_alleles": "contained_mt_allele"}
        )
        bed_temp = bed_temp.annotate(
            locus=bed_temp.contained_mt_allele.locus,
            alleles=bed_temp.contained_mt_allele.alleles,
        )
        bed_temp = bed_temp.key_by(bed_temp.locus, bed_temp.alleles)

        # Annotate back onto the original mt cases where the bed interval falls completely within the variant interval
        mt = mt.annotate_rows(start_and_end_span=bed_temp[mt.locus, mt.alleles].target)

        # Add artifact-prone site filter to any SNP/deletion that starts within, ends within, or completely overlaps an artifact-prone site
        mt = mt.annotate_rows(
            filters=hl.if_else(
                (hl.len(mt.start_overlaps) > 0)
                | (hl.len(mt.end_overlaps) > 0)
                | (hl.is_defined(mt.start_and_end_span)),
                {"artifact_prone_site"},
                {"PASS"},
            )
        )

        mt = mt.drop("region", "start_overlaps", "end_overlaps", "start_and_end_span")

        return mt

    # Equivalent to main function in Broad script.
    def combine_vcfs_worker(
        vcf_path_by_sgid: dict[str, str],
        coverage_mt_path: str,
        artifact_prone_sites_path: str,
        checkpoint_prefix: str,
        output_mt_path: str,
        chunk_size: int,
        minimum_homref_coverage: int,
    ):
        """
        vcf_paths: list of vcf paths for each sample group.
        artifact-prone-sites-path: Path to BED file of artifact-prone sites to flag
            in the FILTER column
        tmp_mt_path: Path to write working copy of output mt
        output_mt_path: Path to export final mt. Base also used to generate .vcf.bgz

        """
        from cpg_utils.hail_batch import init_batch

        init_batch()


        logger.info("Combining VCFs...")
        combined_mt_checkpoint = checkpoint_prefix + 'combined_vcfs.mt'
        combined_mt = join_mitochondria_vcfs_into_mt(
            vcf_path_by_sgid, str(checkpoint_prefix), chunk_size
        )
        combined_mt = combined_mt.checkpoint(combined_mt_checkpoint, overwrite=True)

        logger.info("Removing select sample-level filters...")
        combined_mt = remove_genotype_filters(combined_mt)

        logger.info("Determining homoplasmic reference sites...")
        combined_mt = determine_hom_refs(
            combined_mt, coverage_mt_path, minimum_homref_coverage
        )

        logger.info("Applying artifact_prone_site fiter...")
        combined_mt = apply_mito_artifact_filter(combined_mt, artifact_prone_sites_path)

        logger.info("Writing combined MT and VCF...")
        # Set the file names for output files
        out_vcf = re.sub(r"\.mt$", ".vcf.bgz", output_mt_path)

        combined_mt = combined_mt.checkpoint(output_mt_path, overwrite=True)
        # For the VCF output, join FT values by semicolon
        combined_mt = combined_mt.annotate_entries(
            FT=hl.str(";").join(hl.array(combined_mt.FT))
        )

        META_DICT = {
            "filter": {
                "artifact_prone_site": {
                    "Description": "Variant overlaps an artifact-prone site"
                }
            },
            "format": {
                "DP": {
                    "Description": "Depth of coverage",
                    "Number": "1",
                    "Type": "Integer",
                },
                "FT": {
                    "Description": "Sample-level genotype filters",
                    "Number": ".",
                    "Type": "String",
                },
                "HL": {
                    "Description": "Heteroplasmy level",
                    "Number": "1",
                    "Type": "Float",
                },
                "MQ": {
                    "Description": "Mapping quality",
                    "Number": "1",
                    "Type": "Float",
                },
                "TLOD": {
                    "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
                    "Number": "1",
                    "Type": "Float",
                },
            },
        }

        hl.export_vcf(combined_mt, out_vcf, metadata=META_DICT)

    j.call(
        combine_vcfs_worker,
        vcf_path_by_sgid,
        coverage_mt_path,
        artifact_prone_sites_path,
        str(checkpoint_prefix),
        combined_vcf_mt_path,
        chunk_size,
        minimum_homref_coverage,
    )

    return j
