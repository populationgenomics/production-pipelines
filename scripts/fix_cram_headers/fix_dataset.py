"""
Script to apply `fix_one_header.py` to all post-NAGIM CRAMs
"""

import logging

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, fasta_res_group
from cpg_workflows import get_batch
from cpg_workflows.metamist import get_metamist
from cpg_workflows.targets import Dataset

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


script = (to_path(__file__).parent / 'fix_one_header.py').open().read()


dataset = Dataset(name=get_config()['workflow']['dataset'])
b = get_batch(f'Reheader {dataset}')

for i, entry in enumerate(get_metamist().get_sample_entries(dataset.name)):
    sample = dataset.add_sample(
        id=str(entry['id']),
        external_id=str(entry['external_id']),
        meta=entry.get('meta', {}),
    )
    for seq_type in ['genome', 'exome']:
        cram_path = sample.make_cram_path(sequencing_type=seq_type)
        if not cram_path.exists():
            continue

        out_path = cram_path.path.with_suffix('.REHEADERED.cram')
        logging.info(f'#{i+1} {sample} {cram_path} -> {out_path}')

        j = b.new_job(f'Reheader {cram_path} -> {out_path}')
        j.storage('150G' if seq_type == 'genome' else '50G')

        cmd = f"""\
        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path.path)} {j.out_cram}

        cat <<EOT >> fix_one_header.py
        {script}
        EOT

        samtools reheader {j.out_cram} --in-place \
        --command "fix_one_header.py {fasta_res_group(b)['dict']}"
        """
        j.command(
            command(cmd, monitor_space=True, setup_gcp=True, define_retry_function=True)
        )
        b.write_output(j.out_cram, str(out_path.with_suffix('')))


b.run(wait=False)
