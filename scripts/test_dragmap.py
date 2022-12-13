#!/usr/bin/env python3

from cpg_workflows import get_batch
from cpg_workflows.filetypes import BamPath, CramPath
from cpg_workflows.jobs import align
from cpg_workflows.targets import Sample, Dataset

dataset = Dataset('test')
sample = dataset.add_sample(
    'CHM1_CHM13_2',
    alignment_input_by_seq_type={
        'genome': BamPath(
            'gs://cpg-validation-test-upload/CHM1_CHM13_2.bam',
            index_path='gs://cpg-validation-test-upload/CHM1_CHM13_2.bam.bai',
        )
    },
)

jobs = align.align(
    b=get_batch('Test dragmap'),
    sample=sample,
    output_path=CramPath('gs://cpg-validation-test-tmp/test-dragmap/CHM1_CHM13_2.cram'),
)
get_batch().run(wait=False)
