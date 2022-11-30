"""
Script to apply `fix_one_header.py` to all post-NAGIM CRAMs
"""

import logging

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, fasta_res_group, image_path
from cpg_workflows import get_batch
from cpg_workflows.metamist import get_metamist
from cpg_workflows.targets import Dataset

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

UNMASKED_REF_DICT = (
    'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.dict'
)

script = (to_path(__file__).parent / 'fix_one_header.py').open().read()

dataset = Dataset(name=get_config()['workflow']['dataset'])
b = get_batch(f'Reheader {dataset}')
unmasked_dict = b.read_input(UNMASKED_REF_DICT)

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

        out_path = cram_path.path.parent / 'reheadered' / cram_path.path.name
        logging.info(f'#{i+1} {sample} {cram_path} -> {out_path}')

        j = b.new_job(f'Reheader {cram_path} -> {out_path}')
        j.storage('100G' if seq_type == 'genome' else '30G')
        j.image(image_path('samtools'))

        cmd = f"""\
# Retrying copying to avoid google bandwidth limits
CRAM=$BATCH_TMPDIR/sample.cram
CRAI=$BATCH_TMPDIR/sample.cram.crai
retry_gs_cp {str(cram_path.path)} $CRAM
retry_gs_cp {str(cram_path.path)}.crai $CRAI

cat <<EOT >> fix_one_header.py
{script.replace('`', '')}
EOT

samtools reheader $CRAM --in-place \
--command "python fix_one_header.py {unmasked_dict}"

mv $CRAM {j.out_cram}
mv $CRAM {j.out_crai}
"""
        j.command(
            command(
                cmd,
                monitor_space=True,
                setup_gcp=True,
                define_retry_function=True,
                rm_leading_space=False,
            )
        )
        b.write_output(j.out_cram, str(out_path))
        b.write_output(j.out_crai, str(out_path.with_suffix('.cram.crai')))


b.run(wait=False)
