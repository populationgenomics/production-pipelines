"""
Utilities and resources for running benchmarking
"""
from . import to_path
from .types import FastqPair, CramPath, FastqPairs

BENCHMARK_BUCKET = to_path('gs://cpg-fewgenomes-test/benchmark')
TOY_INPUTS_BUCKET = BENCHMARK_BUCKET / 'inputs/toy'
RESULTS_BUCKET = BENCHMARK_BUCKET / 'outputs'


# 40k reads:
tiny_fq = FastqPairs([
    # This set is 50MB each:
    FastqPair(
        TOY_INPUTS_BUCKET / '2-699835.L001.R1.n40000.fastq.gz',
        TOY_INPUTS_BUCKET / '2-699835.L001.R2.n40000.fastq.gz',
    ),
    FastqPair(
        TOY_INPUTS_BUCKET / '2-699835.L002.R1.n40000.fastq.gz',
        TOY_INPUTS_BUCKET / '2-699835.L002.R2.n40000.fastq.gz',
    ),
])

# ~300k reads:
tiny_cram = CramPath(TOY_INPUTS_BUCKET / 'NA12878-chr21-tiny.cram')

# WGS:
giab_crams = {
    sn: CramPath(f'gs://cpg-reference/validation/giab/cram/{sn}.cram')
    for sn in ['NA12878', 'NA12891', 'NA12892']
}
na12878fq = FastqPairs([
    FastqPair(
        BENCHMARK_BUCKET / 'inputs/NA12878/ERR194147_1.fastq.gz',
        BENCHMARK_BUCKET / 'inputs/NA12878/ERR194147_2.fastq.gz',
    )
])
perth_neuro_fq = FastqPairs([
    FastqPair(
        BENCHMARK_BUCKET
        / 'inputs/PERTHNEURO_FQ/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R1.fastq.gz',
        BENCHMARK_BUCKET
        / 'inputs/PERTHNEURO_FQ/HNFWKCCXY_3_181017_FD07777491_Homo-sapiens__R_170503_GINRAV_DNA_M002_R2.fastq.gz',
    )
])
perth_neuro_cram = CramPath(
    BENCHMARK_BUCKET / 'inputs/PERTHNEURO_CRAM/CPG13045.cram',
)
