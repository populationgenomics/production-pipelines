import collections
import logging

import hail as hl

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build
from cpg_workflows.inputs import get_multicohort
from cpg_workflows.targets import SequencingGroup
from cpg_workflows.utils import can_reuse, exists


def _check_gvcfs(sequencing_groups: list[SequencingGroup]) -> list[SequencingGroup]:
    """
    Making sure each sequencing group has a GVCF
    """
    for sequencing_group in sequencing_groups:
        if not sequencing_group.gvcf:
            if get_config()['workflow'].get('skip_sgs_with_missing_input', False):
                logging.warning(f'Skipping {sequencing_group} which is missing GVCF')
                sequencing_group.active = False
                continue
            else:
                raise ValueError(
                    f'Sequencing group {sequencing_group} is missing GVCF. '
                    f'Use workflow/skip_sgs = [] or '
                    f'workflow/skip_sgs_with_missing_input '
                    f'to control behaviour',
                )

        if get_config()['workflow'].get('check_inputs', True):
            if not exists(sequencing_group.gvcf.path):
                if get_config()['workflow'].get('skip_sgs_with_missing_input', False):
                    logging.warning(f'Skipping {sequencing_group} that is missing GVCF {sequencing_group.gvcf.path}')
                    sequencing_group.active = False
                else:
                    raise ValueError(
                        f'Sequencing group {sequencing_group} is missing GVCF. '
                        f'Use workflow/skip_sgs = [] or '
                        f'workflow/skip_sgs_with_missing_input '
                        f'to control behaviour',
                    )
    return [s for s in sequencing_groups if s.active]


def check_duplicates(iterable):
    """
    Throws error if input list contains repeated items.
    """
    duplicates = [item for item, count in collections.Counter(iterable).items() if count > 1]
    if duplicates:
        raise ValueError(f'Found {len(duplicates)} duplicates: {duplicates}')
    return duplicates


def run(out_vds_path: Path, tmp_prefix: Path, *sequencing_group_ids) -> hl.vds.VariantDataset:
    """
    run VDS combiner, assuming we are on a cluster.
    @param out_vds_path: output path for VDS
    @param tmp_prefix: tmp path for intermediate fields
    @param sequencing_group_ids: optional list of sequencing groups to subset from get_multicohort()
    @return: VDS object
    """
    if can_reuse(out_vds_path):
        return hl.vds.read_vds(str(out_vds_path))

    sequencing_groups = get_multicohort().get_sequencing_groups()
    if sequencing_group_ids:
        sequencing_groups = [s for s in sequencing_groups if s in sequencing_group_ids]

    sequencing_groups = _check_gvcfs(sequencing_groups)

    params = get_config().get('large_cohort', {}).get('combiner', {})

    if intervals := params.get('intervals'):
        if isinstance(intervals, list):
            params['intervals'] = hl.eval(
                [hl.parse_locus_interval(interval, reference_genome=genome_build()) for interval in intervals],
            )
        else:
            params['intervals'] = hl.import_locus_intervals(params['intervals'])
    elif get_config()['workflow']['sequencing_type'] == 'exome':
        params.setdefault('use_exome_default_intervals', True)
    elif get_config()['workflow']['sequencing_type'] == 'genome':
        params.setdefault('use_genome_default_intervals', True)
    else:
        raise ValueError(
            'Either combiner/intervals must be set, or workflow/sequencing_type must be one of: "exome", "genome"',
        )

    sequencing_group_names = [s.id for s in sequencing_groups]
    logging.info(
        f'Combining {len(sequencing_groups)} sequencing groups: '
        f'{", ".join(sequencing_group_names)}, using parameters: '
        f'{params}',
    )

    gvcf_paths = [str(s.gvcf.path) for s in sequencing_groups if s.gvcf]
    if not gvcf_paths:
        raise ValueError('No sequencing groups with GVCFs found')

    logging.info(f'Combining {len(sequencing_group_names)} sequencing groups: {", ".join(sequencing_group_names)}')

    check_duplicates(sequencing_group_names)
    check_duplicates(gvcf_paths)

    combiner = hl.vds.new_combiner(
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sequencing_group_names,
        # Header must be used with gvcf_sample_names, otherwise gvcf_sample_names
        # will be ignored. The first gvcf path works fine as a header because it will
        # be only read until the last line that begins with "#":
        gvcf_external_header=gvcf_paths[0],
        output_path=str(out_vds_path),
        reference_genome='GRCh38',
        temp_path=str(tmp_prefix),
        force=True,
        **params,
    )
    combiner.run()
    return hl.vds.read_vds(str(out_vds_path))
