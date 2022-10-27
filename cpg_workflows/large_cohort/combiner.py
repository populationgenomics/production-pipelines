import collections
import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build
import hail as hl

from cpg_workflows.inputs import get_cohort
from cpg_workflows.utils import can_reuse, exists


def _check_gvcfs():
    """
    Making sure each sample has a GVCF
    """
    for sample in get_cohort().get_samples():
        if not sample.gvcf:
            if get_config()['workflow'].get('skip_samples_with_missing_input', False):
                logging.warning(f'Skipping {sample} which is missing GVCF')
                sample.active = False
                continue
            else:
                raise ValueError(
                    f'Sample {sample} is missing GVCF. '
                    f'Use workflow/skip_samples = [] or '
                    f'workflow/skip_samples_with_missing_input '
                    f'to control behaviour'
                )

        if get_config()['workflow'].get('check_inputs', True):
            if not exists(sample.gvcf.path):
                if get_config()['workflow'].get(
                    'skip_samples_with_missing_input', False
                ):
                    logging.warning(
                        f'Skipping {sample} that is missing GVCF {sample.gvcf.path}'
                    )
                    sample.active = False
                else:
                    raise ValueError(
                        f'Sample {sample} is missing GVCF. '
                        f'Use workflow/skip_samples = [] or '
                        f'workflow/skip_samples_with_missing_input '
                        f'to control behaviour'
                    )


def check_duplicates(iterable):
    """
    Throws error if input list contains repeated items.
    """
    duplicates = [
        item for item, count in collections.Counter(iterable).items() if count > 1
    ]
    if duplicates:
        raise ValueError(f'Found {len(duplicates)} duplicates: {duplicates}')
    return duplicates


def run(out_vds_path: Path, tmp_prefix: Path) -> hl.vds.VariantDataset:
    if can_reuse(out_vds_path):
        return hl.vds.read_vds(str(out_vds_path))

    _check_gvcfs()
    params = get_config().get('combiner', {})

    if intervals := params.get('intervals'):
        if isinstance(intervals, list):
            params['intervals'] = hl.eval(
                [
                    hl.parse_locus_interval(interval, reference_genome=genome_build())
                    for interval in intervals
                ]
            )
        else:
            params['intervals'] = hl.import_locus_intervals(params['intervals'])
    elif get_config()['workflow']['sequencing_type'] == 'exome':
        params.setdefault('use_exome_default_intervals', True)
    elif get_config()['workflow']['sequencing_type'] == 'genome':
        params.setdefault('use_genome_default_intervals', True)
    else:
        raise ValueError(
            'Either combiner/intervals must be set, or workflow/sequencing_type '
            'must be one of: "exome", "genome"'
        )

    logging.info(
        f'Combining {len(get_cohort().get_samples())} samples: '
        f'{", ".join(get_cohort().get_sample_ids())}, using parameters: '
        f'{params}'
    )

    sample_names = get_cohort().get_sample_ids()
    gvcf_paths = [str(s.gvcf.path) for s in get_cohort().get_samples() if s.gvcf]
    logging.info(f'Combining {len(sample_names)} samples: {", ".join(sample_names)}')

    check_duplicates(sample_names)
    check_duplicates(gvcf_paths)

    combiner = hl.vds.new_combiner(
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sample_names,
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
