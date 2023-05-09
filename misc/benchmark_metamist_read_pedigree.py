#!/usr/bin/env python3

"""
Test Metamist's get_pedigree which might get stuck for too many samples, as it uses
a GET request with a limit on the URL string length. This script was used to determine
parameters for the hack in `cpg_workflows.metamist.get_ped_entries()` to get around
this issue.
"""

import logging
import time

from metamist.apis import FamilyApi

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


metamist_proj = 'thousand-genomes'

fapi = FamilyApi()

entries = fapi.get_families(metamist_proj)
family_ids = [entry['id'] for entry in entries]
print(f'families: {family_ids}')
print(f'total families count: {len(family_ids)}')

subset_n = 550
while True:
    start_time = time.time()

    # subset_n = 2**i
    if subset_n >= len(family_ids):
        break

    print(f'Taking first {subset_n} families... ', end='')
    fids = family_ids[:subset_n]
    ped_entries = fapi.get_pedigree(
        internal_family_ids=fids,
        export_type='json',
        project=metamist_proj,
        replace_with_participant_external_ids=True,
    )

    current_time = time.time()
    elapsed_time = current_time - start_time
    print(f'finished in {elapsed_time} seconds')

    subset_n += 10
    time.sleep(0.01)
