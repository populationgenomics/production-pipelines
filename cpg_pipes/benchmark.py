from cpg_pipes.hailbatch import AlignmentInput


BENCHMARK_BUCKET = 'gs://cpg-fewgenomes-test/benchmark'
TOY_INPUTS_BUCKET = f'{BENCHMARK_BUCKET}/inputs/test'
RESULTS_BUCKET = f'{BENCHMARK_BUCKET}/fastqc'


na12878fq = AlignmentInput(
    fqs1=[f'{BENCHMARK_BUCKET}/inputs/NA12878/ERR194147_1.fastq.gz'],
    fqs2=[f'{BENCHMARK_BUCKET}/inputs/NA12878/ERR194147_2.fastq.gz'],
)

tiny_fq = AlignmentInput(
    # This set is 50MB each:
    fqs1=[f'{BENCHMARK_BUCKET}/inputs/toy/2-699835.L001.R1.n40000.fastq.gz'],
    fqs2=[f'{BENCHMARK_BUCKET}/inputs/toy/2-699835.L001.R2.n40000.fastq.gz'],
)

tiny_cram = AlignmentInput(
    bam_or_cram_path=f'{BENCHMARK_BUCKET}/inputs/toy/NA12878-chr21-tiny.cram',
)

giab_cram_inputs = {
    sn: AlignmentInput(
        bam_or_cram_path=f'gs://cpg-reference/validation/giab/cram/{sn}.cram'
    )
    for sn in ['NA12878', 'NA12891', 'NA12892']
}
