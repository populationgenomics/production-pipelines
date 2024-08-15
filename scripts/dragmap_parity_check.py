#! /usr/bin/env python3
import json
import subprocess
from os.path import basename

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, get_batch, init_batch
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

    out_path = f'gs://cpg-{project}/dragmap_parity/{"nagim_" if nagim else "new_"}vds.vds'
    # make Combiner objects
    combiner = hl.vds.new_combiner(
        output_path=out_path,
        temp_path=f'gs://cpg-{project}-tmp/dragmap_parity/',
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sample_names,
        gvcf_external_header=external_header,
        reference_genome='GRCh38',
        use_genome_default_intervals=True,
    )
    return combiner, out_path


def create_vds(project: str) -> hl.vds.VariantDataset:
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
    nagim_combiner, nagim_vds_path = create_combiner(
        project=project,
        gvcf_paths=checked_nagim_gvcf_paths,
        sample_names=expids,
        external_header=checked_nagim_gvcf_paths[0],
        nagim=True,
    )
    new_combiner, new_vds_path = create_combiner(
        project=project,
        gvcf_paths=checked_new_gvcf_paths,
        sample_names=expids,
        external_header=checked_new_gvcf_paths[0],
    )

    # run the combiners
    nagim_combiner.run()
    new_combiner.run()

    return hl.vds.read_vds(nagim_vds_path), hl.vds.read_vds(new_vds_path)


@click.command()
@click.option('--project', required=True)
def main(project: str):

    nagim_vds, new_vds = create_vds(project)

    # prepare vds' for comparison
    # As per documentation, hl.methods.concordance() requires the dataset to contain no multiallelic variants.
    # as well as the entry field to be 'GT', also expects MatrixTable, not a VariantDataset
    nagim_vds = hl.vds.split_multi(nagim_vds, filter_changed_loci=True)
    new_vds = hl.vds.split_multi(new_vds, filter_changed_loci=True)
    nagim_mt: hl.MatrixTable = nagim_vds.variant_data
    new_mt: hl.MatrixTable = new_vds.variant_data
    nagim_mt = nagim_mt.annotate_entries(GT=hl.vds.lgt_to_gt(nagim_mt.LGT, nagim_mt.LA))
    new_mt = new_mt.annotate_entries(GT=hl.vds.lgt_to_gt(new_mt.LGT, new_mt.LA))
    # checkpoint the MatrixTables
    nagim_mt = nagim_mt.checkpoint(dataset_path('/dragmap_parity/nagim_mt.mt', 'tmp'))
    new_mt = new_mt.checkpoint(dataset_path('/dragmap_parity/new_mt.mt', 'tmp'))

    # compare the two VDS'
    global_conc, cols_conc, rows_conc = hl.concordance(nagim_mt, new_mt)

    # write the concordance results
    with open(dataset_path('/dragmap_parity/global_concordance.json'), 'w') as f:
        json.dump(global_conc, f)
    cols_conc.checkpoint(dataset_path('/dragmap_parity/cols_concordance.ht'))
    rows_conc.checkpoint(dataset_path('/dragmap_parity/rows_concordance.ht'))


if __name__ == '__main__':
    init_batch()
    main()
