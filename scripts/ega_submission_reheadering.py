import logging
from typing import TYPE_CHECKING

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch, init_batch

if TYPE_CHECKING:
    from hailtop.batch.job import Job

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


@click.command()
@click.option(
    '--tenk10k-sgid-mapping-path',
    type=str,
    required=True,
    help='Path to CSV mapping CPG sequencing group IDs to EGA IDs.',
)
@click.option(
    '--cpg-id-cram-mapping-path',
    type=str,
    required=True,
    help='Path to CSV mapping CPG sequencing group IDs to CRAM file paths.',
)
@click.option('--output-dir', type=str, required=True, help='Directory to write reheadered CRAM files.')
@click.option('--dry-run', is_flag=True, help='If set, the batch will not be run.')
def main(
    tenk10k_sgid_mapping_path: str,
    cpg_id_cram_mapping_path: str,
    output_dir: str,
    dry_run: bool,
):
    init_batch()
    b = get_batch()
    ref_fasta = fasta_res_group(b)

    tenk10k_sgid_mapping = pd.read_csv(tenk10k_sgid_mapping_path)
    cpg_id_cram_mapping = pd.read_csv(cpg_id_cram_mapping_path)

    cpg_id_cram_map_dict = cpg_id_cram_mapping.set_index('sequencing_group_id')['cram_filename'].to_dict()
    tenk10k_sgid_map_dict = tenk10k_sgid_mapping.set_index('cpg_sequencing_group_id')['tenk_id'].to_dict()

    jobs: list[Job] = []
    for cpg_id, ega_id in tenk10k_sgid_map_dict.items():
        if cpg_id not in cpg_id_cram_map_dict:
            logger.warning(f'No CRAM found for CPG ID {cpg_id}, skipping...')
            continue
        incram = cpg_id_cram_map_dict[cpg_id]
        out_path_base = str(to_path(output_dir) / f'{ega_id}')
        input_cram_reads = b.read_input_group(
            **{
                'cram': str(incram),
                'cram.crai': str(f'{incram}.crai'),
            },
        )
        j = b.new_bash_job(name=f'reheader_{incram}_to_{out_path_base}.cram')
        j.image(image_path('samtools'))
        j.storage('200Gi')
        j.memory('8Gi')

        j.declare_resource_group(
            output_cram={
                'cram': '{root}.cram',
                'cram.crai': '{root}.cram.crai',
            },
        )
        j.command(
            f"""
            set -e

            # Explanation of sed commands:
            # 1. /^@RG/   -> On Read Group lines: Replace specific CPG ID with EGA ID
            # 2. /^[^@]/  -> On Body lines (not starting with @): Replace specifically the 'RG:Z:' tag
            # 3. /^@PG/   -> On Program History lines: Redact ANY CPG-like ID (e.g., CPGXXXXXX) with 'REDACTED'

            samtools view -h -T "{ref_fasta.base}" "{input_cram_reads.cram}" \\
            | sed -e '/^@RG/ s/{cpg_id}/{ega_id}/g' \\
                  -e '/^[^@]/ s/RG:Z:{cpg_id}/RG:Z:{ega_id}/g' \\
                  -e '/^@PG/ s/CPG[0-9]*/REDACTED/g' \\
            | samtools view -C -T "{ref_fasta.base}" --no-PG -o "{j.output_cram.cram}"

            samtools index "{j.output_cram.cram}" "{j.output_cram['cram.crai']}"
            """,
        )

        b.write_output(j.output_cram, out_path_base)
        jobs.append(j)

    b.run(wait=False, dry_run=dry_run)


if __name__ == '__main__':
    main()
