import subprocess
from os.path import basename

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path
from metamist.graphql import gql, query

ACTIVE_INACTIVE_QUERY = gql(
    """
    query getActiveInactiveSGs($project: String!) {
        project(name: $project) {
            id
            samples {
                externalId
                id
                inactive_sgs: sequencingGroups(activeOnly: {eq: false}) {
                    id
                }
                active_sgs: sequencingGroups(activeOnly: {eq: true}) {
                    id
                }
            }
        }
    }
    """,
)

config = get_config()


def get_dict_of_gvcf_directories(project: str, nagim: bool = False) -> dict[str, str]:
    sgid_gvcf_path_map = {}
    prefix = f'gs://cpg-{project}/gvcf/{"nagim" if nagim else ""}'

    # Get GCP file paths
    cmd = ['gsutil', 'ls', prefix]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        gvcf_paths: list[str] = result.stdout.splitlines()
        for gvcf_path in gvcf_paths:
            if gvcf_path.endswith('.gz'):
                sgid = basename(gvcf_path).split('.')[0]
                sgid_gvcf_path_map[sgid] = gvcf_path
    return sgid_gvcf_path_map


def create_combiner(
    project: str,
    gvcf_paths: list[str],
    sample_names: list[str],
    external_header: str,
    nagim: bool = False,
) -> hl.vds.combiner.VariantDatasetCombiner:
    # make Combiner objects
    combiner = hl.vds.new_combiner(
        output_path=f'gs://cpg-{project}/dragmap_parity/{"nagim_" if nagim else "new_"}vds.vds',
        temp_path=f'gs://cpg-{project}-tmp/dragmap_parity/',
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sample_names,
        gvcf_external_header=external_header,
        reference_genome='GRCh38',
        use_genome_default_intervals=True,
    )
    return combiner


def parity_check(project: str):
    active_inactive_response = query(ACTIVE_INACTIVE_QUERY, variables={'project': project})

    active_inactive_sg_map = {}
    for sample in active_inactive_response['project']['samples']:
        if len(sample['inactive_sgs']) >= 1:
            assert len(sample['active_sgs']) == len(sample['inactive_sgs']) == 1
            active_inactive_sg_map[sample['externalId']] = {
                'active': sample['active_sgs'][0]['id'],
                'inactive': sample['inactive_sgs'][0]['id'],
            }
    nagim_gvcf_paths_dict: dict[str, str] = get_dict_of_gvcf_directories(project, nagim=True)
    gvcf_paths_dict: dict[str, str] = get_dict_of_gvcf_directories(project)

    nagim_sgids = set(nagim_gvcf_paths_dict.keys())
    new_sgids = set(gvcf_paths_dict.keys())

    checked_nagim_gvcf_paths = []
    checked_new_gvcf_paths = []
    expids = []
    for expid, sgids in active_inactive_sg_map.items():
        inactive_sgids = sgids['inactive']
        active_sgid = sgids['active']
        # check inactive sgid has nagim gvcf
        if inactive_sgids in nagim_sgids and active_sgid in new_sgids:
            checked_nagim_gvcf_paths.append(nagim_gvcf_paths_dict[inactive_sgids])
            checked_new_gvcf_paths.append(gvcf_paths_dict[active_sgid])
            expids.append(expid)

    # make Combiner objects
    nagim_combiner: hl.vds.combiner.VariantDatasetCombiner = create_combiner(
        project=project,
        gvcf_paths=checked_nagim_gvcf_paths,
        sample_names=expids,
        external_header=checked_nagim_gvcf_paths[0],
        nagim=True,
    )
    new_combiner: hl.vds.combiner.VariantDatasetCombiner = create_combiner(
        project=project,
        gvcf_paths=checked_new_gvcf_paths,
        sample_names=expids,
        external_header=checked_new_gvcf_paths[0],
    )

    # run the combiners
    nagim_combiner.run()
    new_combiner.run()


@click.option('--project', required=True)
def main(project: str):
    b = get_batch('Dragmap parity check')

    j = b.new_python_job(name='Dragmap Parity Check')
    j.image(config['workflow']['driver_image'])
    j.call(parity_check, project)

    b.run(wait=False)


if __name__ == '__main__':
    main()
