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
    tsv_string = ""
    for sid, path in base_level_coverage_by_sid.items():
        tsv_string += f'{sid}\\t{path}\\t{sid}\\n'

    # script needs to write to a .ht path
    j.declare_resource_group(
        outfile={'ht': '{root}.ht'}
    )
    cmd = f"""
        # build input file:
        echo "{tsv_string}" > input.tsv
        cp input.tsv {j.tsv}

        ls -l

        ls -l /

        # Run query job
        # python {annotate_coverage_script.__file__}
        # python cpg_workflows/mito_pipeline_scripts/annotate_coverage.py \
        #     --input-tsv input.tsv \
        #     --output-ht {j.outfile.ht} \
        #     --temp-dir $BATCH_TMPDIR/mt
        """

    j.command(command(cmd, setup_gcp=True, monitor_space=True))
    # b.write_output(j.outfile, str(coverage_ht.with_suffix('')))
    b.write_output(j.tsv, str(coverage_ht.with_suffix('')) + '.temp.tsv')

    return j
