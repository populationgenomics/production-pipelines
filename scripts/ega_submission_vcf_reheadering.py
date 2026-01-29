from operator import index
from typing import TYPE_CHECKING

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch, init_batch
from cpg_workflows.jobs import vcf

if TYPE_CHECKING:
    from hailtop.batch.job import Job


@click.command()
@click.option(
    '--tenk10k-sgid-mapping-path',
    type=str,
    required=True,
    help='Path to CSV mapping CPG sequencing group IDs to EGA IDs.',
)
@click.option('--tr-vcfs', is_flag=True, help='If set, cut "CPG" prefix from CPG IDs in tenk10k mapping.')
@click.option(
    '--tr-vcf-path-csv',
    type=str,
    help='Path to CSV containing VCF paths that are from Hope\'s TR VCF cohort.',
)
@click.option(
    '--com-rare-vcf-path-csv',
    type=str,
    help='Path to CSV containing VCF paths that are from Anna\'s common/rare VCFs.',
)
@click.option('--output-dir', type=str, required=True, help='Directory to write reheadered CRAM files.')
@click.option('--dry-run', is_flag=True, help='If set, the batch will not be run.')
def main(
    tenk10k_sgid_mapping_path: str,
    tr_vcfs: bool,
    tr_vcf_path_csv: str,
    com_rare_vcf_path_csv: str,
    output_dir: str,
    dry_run: bool,
):
    if tr_vcf_path_csv and not tr_vcfs:
        raise ValueError('If --tr-vcf-path-csv is provided, --tr-vcfs must also be set.')

    init_batch()
    b = get_batch()

    # read mapping files
    tenk10k_sgid_mapping = pd.read_csv(tenk10k_sgid_mapping_path)

    tenk10k_sgid_map_dict = tenk10k_sgid_mapping.set_index('cpg_sequencing_group_id')['tenk_id'].to_dict()

    if tr_vcfs:
        index_ext = 'tbi'
        index_flag = '-t'
        path_df = pd.read_csv(tr_vcf_path_csv)
        # cut 'CPG' prefix from CPG IDs in tenk10k_sgid_map_dict keys if Hope's VCFs are used
        # map cpg IDs without 'CPG' prefix to tenk10k IDs
        tenk10k_sgid_map_dict = {k[3:]: v for k, v in tenk10k_sgid_map_dict.items()}
    else:
        index_ext = 'csi'
        index_flag = '-c'
        path_df = pd.read_csv(com_rare_vcf_path_csv)

    mapping_lines = []
    for cpg_id, tenk_id in tenk10k_sgid_map_dict.items():
        # Handle the CPG prefix logic here if needed, or assume map dict is already clean
        mapping_lines.append(f"{cpg_id} {tenk_id}")

    # Create the block of text for the bash file
    mapping_block = '\n'.join(mapping_lines)

    jobs: list[Job] = []
    for i, vcf_path in enumerate(path_df['vcf_path'].tolist()):

        out_path_base = str(to_path(output_dir) / f'{to_path(vcf_path).name}_reheadered')

        # Define Input
        input_vcf = b.read_input_group(
            **{'vcf': str(vcf_path), 'vcf.csi': str(f'{vcf_path}.{index_ext}')},
        )

        # Create ONE Job
        j = b.new_bash_job(name=f'process_{vcf_path}')
        j.image(image_path('bcftools'))
        j.storage('100Gi')
        j.memory('4Gi')

        j.declare_resource_group(
            # can't use f-strings here due to hail batch parsing
            output_vcf={'vcf': '{root}.vcf.gz', 'vcf.csi': '{root}.vcf.gz.' + index_ext},
        )

        j.command(
            f"""
            set -e

            # 1. Create the master mapping file (OldName NewName)
            cat <<EOF > rename_map.txt
{mapping_block}
EOF

            # 2. Extract just the OldNames to create the subset list
            awk '{{print $1}}' rename_map.txt > subset_list.txt

            # 3. The Pipeline: Subset -> Filter -> Reheader -> Write
            # -Ou : Output uncompressed BCF (fastest for piping)
            # --threads 4 : Speed up the compression at the end

            bcftools view -S subset_list.txt --force-samples "{input_vcf.vcf}" -Ou \\
            | bcftools view -c 1 -a -Ou \\
            | bcftools reheader -s rename_map.txt --output "{j.output_vcf.vcf}"

            # 4. Index the final result
            bcftools index {index_flag} "{j.output_vcf.vcf}"
            """,
        )

        b.write_output(j.output_vcf, out_path_base)

        jobs.append(j)

    b.run(wait=False, dry_run=dry_run)


if __name__ == '__main__':
    main()
