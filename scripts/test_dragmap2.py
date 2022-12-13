#!/usr/bin/env python3

from cpg_workflows import get_batch
from cpg_workflows.filetypes import CramPath, FastqPair, FastqPairs
from cpg_workflows.jobs import align
from cpg_workflows.targets import Dataset

dataset = Dataset('test')
sample = dataset.add_sample(
    'FROM_FASTQ',
    alignment_input_by_seq_type={
        'genome': FastqPairs(
            [
                FastqPair(
                    'gs://cpg-validation-test-upload/HCMVGDSX3_1_220405_FD07777372_Homo-sapiens_TCCGCCAATT-CAGCACGGAG_R_220405_CNTROL_DNA_M001_R1.fastq.gz',
                    'gs://cpg-validation-test-upload/HCMVGDSX3_1_220405_FD07777372_Homo-sapiens_TCCGCCAATT-CAGCACGGAG_R_220405_CNTROL_DNA_M001_R2.fastq.gz',
                )
            ]
        )
    },
)

jobs = align.align(
    b=get_batch('Test dragmap from fq'),
    sample=sample,
    output_path=CramPath('gs://cpg-validation-test-tmp/test-dragmap/from_fq.cram'),
)
get_batch().run(wait=False)
