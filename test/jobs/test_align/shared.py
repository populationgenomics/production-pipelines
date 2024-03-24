from typing import Literal

from hailtop.batch.job import Job

from ...factories.config import PipelineConfig, WorkflowConfig


def default_config() -> PipelineConfig:
    return PipelineConfig(
        workflow=WorkflowConfig(
            dataset='align-test',
            access_level='test',
            sequencing_type='genome',
            check_inputs=False,
        ),
        images={
            'dragmap': 'dragmap:latest',
            'bwa': 'bwa:latest',
            'bwamem2': 'bwamem2:latest',
            'samtools': 'samtools:latest',
            'picard': 'picard:latest',
        },
        references={
            'broad': {
                'ref_fasta': 'broad_reference.fa',
                'dragmap_prefix': 'a-cpg-bucket/dragen_reference/',
            },
        },
        other={'resource_overrides': {}},
    )


def select_jobs(jobs: list[Job], jtype: Literal['align', 'merge', 'markdup']) -> list[Job]:
    if jtype == 'align':
        return [j for j in jobs if 'Align' in str(j.name)]
    elif jtype == 'merge':
        return [j for j in jobs if 'Merge' in str(j.name)]
    elif jtype == 'markdup':
        return [j for j in jobs if 'MarkDuplicates' in str(j.name)]
    else:
        raise ValueError(f'Unknown job type: {jtype}')
