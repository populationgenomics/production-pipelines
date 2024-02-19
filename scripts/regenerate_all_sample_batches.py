#!/usr/bin/env python3


"""
Aim is to take all previously batched samples and re-batch them to
get the best possible groupings of samples for GATK-SV

To do this we join every QC table from Evidence QC which has run so far,
filter to retain only active SG IDs, then re-batch cleanly.
"""

import json
import pandas as pd
from metamist.graphql import gql, query

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.jobs.sample_batching import batch_sgs


