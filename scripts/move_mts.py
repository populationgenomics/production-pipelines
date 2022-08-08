"""
Copy dataset mts from seqr bucket to dataset buckets.
"""

from cpg_pipes.hailquery import init_batch
import hail as hl

init_batch()

for path in [
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-acute-care.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-circa.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-heartkids.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-hereditary-neuro.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-mito-disease.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-ohmr3-mendelian.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-ohmr4-epilepsy.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-perth-neuro.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-ravenscroft-arch.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-ravenscroft-rdstudy.mt',
    'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-validation.mt',
]:
    mt = hl.read_matrix_table(path)
    dataset = path.replace(
        'gs://cpg-seqr-main/mt/e51f4fb948f27a4130f4a56b32fd1ca8e7c0ad_867-', ''
    ).replace('.mt', '')
    new_path = path.replace('gs://cpg-seqr-main', f'gs://cpg-{dataset}-main')
    print(f'{dataset}: moving {path} to {new_path}')
    mt.write(new_path)
