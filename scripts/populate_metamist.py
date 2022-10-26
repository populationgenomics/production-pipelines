"""
Script to populate Metamist Analysis entries with status=completed. 

Useful to back-populate entries for old workflow runs. For new runs, TOML entry
`workflow/status_reporter = metamist` should make sure the workflow does populate 
all the entries automatically. But if it wasn't enabled during the workflow run, 
this script would back-populate for an already finished workflow.
"""
import json
import logging
import os
from collections import defaultdict

import click

from cpg_utils.config import get_config, set_config_paths
from inputs import get_cohort
from status import MetamistStatusReporter
from utils import exists

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


COMMANDS = ['analyses', 'qc', 'joint-calling', 'es-index']


@click.command()
@click.argument('command', type=click.Choice(COMMANDS))
@click.option('--config', 'config_paths', multiple=True)
def main(command: str, config_paths: list[str]):
    if _env_var := os.environ.get('CPG_CONFIG_PATH'):
        config_paths += _env_var.split(',') + list(config_paths)
    set_config_paths(list(config_paths))

    sequencing_type = get_config()['workflow']['sequencing_type']
    cohort = get_cohort()
    status = MetamistStatusReporter()

    if command == 'analyses':
        from sample_metadata.apis import AnalysisApi
        from sample_metadata.models import AnalysisQueryModel, AnalysisStatus

        analyses = AnalysisApi().query_analyses(
            AnalysisQueryModel(
                projects=[d.name for d in cohort.get_datasets()],
                status=AnalysisStatus('completed'),
            )
        )
        existing_paths = set(a['output'] for a in analyses)

        for i, sample in enumerate(cohort.get_samples()):
            if path := sample.make_cram_path().path:
                if str(path) in existing_paths:
                    continue

                if not path.exists():
                    continue

                print(f'#{i+1} {sample} {path}')
                status.create_analysis(
                    str(path),
                    analysis_type='cram',
                    analysis_status='completed',
                    target=sample,
                    meta=sample.get_job_attrs()
                    | dict(
                        size=path.stat().st_size,
                        sequencing_type=sequencing_type,
                    ),
                    project_name=sample.dataset.name,
                )
            if (path := sample.make_gvcf_path().path).exists():
                print(f'#{i+1} {sample} {path}')
                if str(path) in existing_paths:
                    continue
                status.create_analysis(
                    str(path),
                    analysis_type='gvcf',
                    analysis_status='completed',
                    target=sample,
                    meta=sample.get_job_attrs()
                    | dict(
                        size=path.stat().st_size,
                        sequencing_type=sequencing_type,
                    ),
                    project_name=sample.dataset.name,
                )

    if command == 'qc':
        """
        Add Analysis entries of type="qc" from a MultiQC JSON data.
        """
        multiqc_json_path = (
            cohort.analysis_dataset.prefix() / 'qc' / 'multiqc_data.json'
        )
        assert multiqc_json_path.exists()

        with multiqc_json_path.open() as f:
            data = json.load(f)

        metrics_by_sample: defaultdict[str, dict] = defaultdict()
        for sample_d in data['report_general_stats_data']:
            for sid, metrics_d in sample_d.items():
                if sid not in metrics_by_sample:
                    metrics_by_sample[sid] = metrics_d
                metrics_by_sample[sid] |= metrics_d

        for i, sample in enumerate(cohort.get_samples()):
            print(f'#{i+1} {sample}')
            if sample.rich_id not in metrics_by_sample:
                print(f'{sample.rich_id} not found in MultiQC, skipping')
                continue
            metrics_d = metrics_by_sample[sample.rich_id]
            status.create_analysis(
                str(multiqc_json_path),
                analysis_type='qc',
                analysis_status='completed',
                target=sample,
                meta=sample.get_job_attrs()
                | dict(
                    sequencing_type=sequencing_type,
                    metrics=metrics_d,
                ),
                project_name=sample.dataset.name,
            )

    if command == 'joint-calling':
        h = cohort.alignment_inputs_hash()
        path = cohort.analysis_dataset.prefix() / 'mt' / f'{h}.mt'
        if exists(path):
            status.create_analysis(
                str(path),
                analysis_type='joint-calling',
                analysis_status='completed',
                target=cohort,
                meta=cohort.get_job_attrs()
                | dict(
                    sequencing_type=sequencing_type,
                ),
                project_name='seqr',
            )

    if command == 'es-index':
        names = []
        if sequencing_type == 'genome':
            names = [
                'acute-care-genome-2022_0815_1644_xkhvx',
                'validation-genome-2022_0810_2358_474tt',
                'ravenscroft-arch-genome-2022_0618_1137_4qfyn',
                'circa-genome-2022_0618_1137_4qfyn',
                'ohmr3-mendelian-genome-2022_0618_1137_4qfyn',
                'mito-disease-genome-2022_0618_1137_4qfyn',
                'perth-neuro-genome-2022_0618_1137_4qfyn',
                'ohmr4-epilepsy-genome-2022_0618_1137_4qfyn',
                'hereditary-neuro-genome-2022_0618_1137_4qfyn',
                'ravenscroft-rdstudy-genome-2022_0618_1137_4qfyn',
                'heartkids-genome-2022_0618_1137_4qfyn',
            ]
        if sequencing_type == 'exome':
            names = [
                'acute-care-exome-2022_0815_2246_h2pb9',
                'mito-disease-exome-2022_0815_2246_h2pb9',
                'hereditary-neuro-exome-2022_0815_2246_h2pb9',
                'kidgen-exome-2022_0815_2246_h2pb9',
            ]
        for name in names:
            ds_name = name.split(f'-{sequencing_type}-')[0]
            print(f'Adding {ds_name}')
            dataset = cohort.create_dataset(ds_name)
            status.create_analysis(
                str(name),
                analysis_type='es-index',
                analysis_status='completed',
                target=dataset,
                meta=dataset.get_job_attrs()
                | dict(
                    sequencing_type=sequencing_type,
                ),
                project_name=ds_name,
            )


if __name__ == '__main__':
    main()
