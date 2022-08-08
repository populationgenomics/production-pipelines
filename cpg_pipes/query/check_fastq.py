import pysam


def _main(fq_path):
    with pysam.FastxFile(fq_path) as fh:
        for i, entry in enumerate(fh):
            seq_len = len(entry.sequence)
            qua_len = len(entry.quality)
            if seq_len != qua_len:
                print(
                    f'Entry {i} has mismatching sizes of seq ({seq_len}) '
                    f'and qual ({qua_len}):'
                )
                print(entry)
                print(entry.sequence)
                print(entry.quality)
            if i % 10_000_000 == 0:
                print(f'Line {i}...')
    print(f'Total reads iterated: {i}')
