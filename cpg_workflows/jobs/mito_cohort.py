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

from cpg_workflows.mito_pipeline_scripts import annotate_coverage as annotate_coverage_script


def annotate_coverage(
    b,
    base_level_coverage_by_sid: dict[str, Path],
    coverage_ht: Path,
    # haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Combine individual mitochondria coverage files and outputs a hail table with coverage annotations
    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_annotate_coverage', job_attrs)
    # j.image(image_path('haplocheckcli'))
    j.image(get_config()['workflow']['driver_image'])

    res = STANDARD.request_resources(ncpu=8)
    res.set_to_job(j)

    # generate tsv string to use as input file
    # required format: "tsv of participant_id, base_level_coverage_metrics, sample"
    tsv_string = "#participant_id\\tpath\\tsample_id"
    for sid, path in base_level_coverage_by_sid.items():
        tsv_string += f'{sid}\\t{path}\\t{sid}\\n'

    # script needs to write to a .ht path
    j.declare_resource_group(
        outfile={'ht': '{root}.ht'}
    )
    cmd = f"""
        # build input file:
        printf "{tsv_string}" > input.tsv
        cp input.tsv {j.tsv}

        mkdir -p $BATCH_TMPDIR/mito/

        # Run query job
        # python {annotate_coverage_script.__file__}
        python cpg_workflows/mito_pipeline_scripts/annotate_coverage.py \
            --input-tsv input.tsv \
            --output-ht {j.outfile.ht} \
            --temp-dir $BATCH_TMPDIR/mito/
        """

    j.command(command(cmd, setup_gcp=True, monitor_space=True))
    b.write_output(j.outfile, str(coverage_ht.with_suffix('')))
    b.write_output(j.tsv, 'gs://cpg-acute-care-test/mito/input.temp.tsv')

    return j


def annotate_coverage2(
    b,
    base_level_coverage_by_sid: dict[str, Path],
    coverage_ht: Path,
    # haplocheck_output: Path | None,
    job_attrs: dict | None = None,
) -> Job:
    """
    Try using python job....
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
    for sid, path in base_level_coverage_by_sid.items():
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
    import os
    import re
    import sys

    import hail as hl

    from os.path import dirname
    from hail.utils.java import info

    logging.basicConfig(
        format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    logger = logging.getLogger("Annotate coverage")
    logger.setLevel(logging.INFO)

    def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int) -> hl.MatrixTable:
        """
        Hierarchically join together MatrixTables in the provided list.

        :param mts: List of MatrixTables to join together
        :param temp_dir: Path to temporary directory for intermediate results
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

                merged.show(10)
                logger.info(f"merged.count(): {merged.count()}")
                # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
                merged = merged.annotate_globals(
                    __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
                )
                merged.show(10)
                next_stage.append(
                    merged.checkpoint(
                        os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
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

    def annotate_coverage_worker(inputs: list, output_ht: str , temp_dir: str = "", chunk_size: int = 100, overwrite=True):

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
        out_dir = dirname(output_ht)

        cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size)
        n_samples = cov_mt.count_cols()

        logger.info("Adding coverage annotations...")
        # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
        cov_mt = cov_mt.annotate_rows(
            locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
            mean=hl.float(hl.agg.mean(cov_mt.coverage)),
            median=hl.median(hl.agg.collect(cov_mt.coverage)),
            over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
            over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
        )
        cov_mt.show()

        cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

        output_mt = re.sub(r"\.ht$", ".mt", output_ht)
        output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
        output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

        logger.info("Writing sample level coverage...")
        sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
        sample_mt.coverage.export(output_samples)

        logger.info("Writing coverage mt and ht...")
        cov_mt.write(output_mt, overwrite=overwrite)
        cov_ht = cov_mt.rows()
        cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
        cov_ht.export(output_tsv)

    ##################
    ##################

    annotate_coverage_out = j.call(
        annotate_coverage_worker, inputs, j.outfile
    )
    j.outfile.add_extension('.ht')

    return j, annotate_coverage_out
