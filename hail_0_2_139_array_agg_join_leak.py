#!/usr/bin/env python3
"""Hail 0.2.139 regression: array_agg over an entry array field, on a re-keyed MatrixTable, joined
with a second row table from the same MatrixTable, allocates off-heap (RegionPool) memory
geometrically inside a single task.

    python hail_0_2_139_array_agg_join_leak.py [n_rows] [leak|no_rekey|no_join|entry_scalar|col_pred]

On Hail 0.2.139 `leak` with 5000 rows reaches 16 GB of RegionPool allocation within about a
second (64 KB regions doubling: 2M, 8M, 32M, ..., 16G) and is OOM-killed on smaller machines.
On Hail 0.2.138 the same script peaks at 1 MB. Each control removes one ingredient and stays
at about 1 MB on 0.2.139.
"""
import os, sys, time
import hail as hl

n_rows = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
mode = sys.argv[2] if len(sys.argv) > 2 else 'leak'
os.environ['PYSPARK_SUBMIT_ARGS'] = '--driver-memory 2g pyspark-shell'
hl.init(master='local[1]', quiet=True, backend='spark', log=f'hail_{mode}_{n_rows}.log')

mt = hl.utils.range_matrix_table(n_rows=n_rows, n_cols=15, n_partitions=1)
if mode != 'no_rekey':
    mt = mt.key_rows_by(k=mt.row_idx + 1)                        # 1. re-key the rows
mt = mt.annotate_entries(LA=hl.range(0, 3), x=mt.col_idx)         # an entry array field (like a VDS's LA)

if mode == 'entry_scalar':
    pred = lambda ai: hl.agg.count_where(mt.x >= ai)              # entry scalar instead of entry array
elif mode == 'col_pred':
    pred = lambda ai: hl.agg.count_where(mt.col_idx >= ai)        # column field, no entry field
else:
    pred = lambda ai: hl.agg.count_where(mt.LA.contains(ai))      # 2. array_agg reading an entry array
ht = mt.annotate_rows(arr=hl.agg.array_agg(pred, hl.range(1, 3))).rows()

if mode != 'no_join':                                             # 3. join with another row table from the same MT
    other = mt.annotate_rows(dp=hl.agg.sum(mt.col_idx)).rows()
    ht = ht.select(arr=ht.arr, dp=other[ht.key].dp)

t0 = time.time()
print(f'{mode} rows={n_rows}:', ht.head(1).collect()[0], f'in {time.time() - t0:.1f}s', flush=True)
