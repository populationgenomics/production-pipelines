import collections
import logging

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build
from cpg_utils.workflows.inputs import get_cohort
from cpg_utils.workflows.utils import can_reuse, exists
import hail as hl

logger = logging.getLogger(__file__)


def _check_gvcfs():
    """
    Making sure each sample has a GVCF
    """
    for sample in get_cohort().get_samples():
        if not sample.gvcf:
            if get_config()['workflow'].get('skip_samples_with_missing_input', False):
                logger.warning(f'Skipping {sample} which is missing GVCF')
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
                    logger.warning(
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

    logger.info(
        f'Combining {len(get_cohort().get_samples())} samples: '
        f'{", ".join(get_cohort().get_sample_ids())}, using parameters: '
        f'{params}'
    )

    sample_names = get_cohort().get_sample_ids()
    gvcf_paths = [str(s.gvcf.path) for s in get_cohort().get_samples()]
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


# def queue_combiner(out_vds_path: Path) -> Job | None:
#     """
#     Add VCF combiner jobs, produce VDC.
#     """
#     if can_reuse(out_vds_path):
#         return None
#
#     # Combiner takes advantage of autoscaling cluster policies
#     # to reduce costs for the work that uses only the driver machine:
#     # https://hail.is/docs/0.2/experimental/vcf_combiner.html#pain-points
#     # To add a 50-worker policy for a project "prophecy-339301":
#     # ```
#     # gcloud dataproc autoscaling-policies import vcf-combiner-50 \
#     # --source=combiner-autoscaling-policy-50.yaml --region=australia-southeast1 \
#     # --project prophecy-339301
#     # ```
#     scatter_count = get_config()['workflow'].get('scatter_count', 50)
#     if scatter_count > 100:
#         autoscaling_workers = '200'
#     elif scatter_count > 50:
#         autoscaling_workers = '100'
#     else:
#         autoscaling_workers = '50'
#
#     for sample in get_cohort().get_samples():
#         if not sample.gvcf:
#             if get_config()['workflow'].get('skip_samples_with_missing_input', False):
#                 logger.warning(f'Skipping {sample} which is missing GVCF')
#                 sample.active = False
#                 continue
#             else:
#                 raise ValueError(
#                     f'Sample {sample} is missing GVCF. '
#                     f'Use workflow/skip_samples = [] or '
#                     f'workflow/skip_samples_with_missing_input '
#                     f'to control behaviour'
#                 )
#
#         gvcf_path = sample.gvcf.path
#         if get_config()['workflow'].get('check_inputs', True):
#             if not exists(gvcf_path):
#                 if get_config()['workflow'].get(
#                     'skip_samples_with_missing_input', False
#                 ):
#                     logger.warning(
#                         f'Skipping {sample} that is missing GVCF {gvcf_path}'
#                     )
#                     sample.active = False
#                 else:
#                     raise ValueError(
#                         f'Sample {sample} is missing GVCF. '
#                         f'Use workflow/skip_samples = [] or '
#                         f'workflow/skip_samples_with_missing_input '
#                         f'to control behaviour'
#                     )
#
#     logger.info(
#         f'Combining {len(get_cohort().get_samples())} samples: '
#         f'{", ".join(get_cohort().get_sample_ids())}'
#     )
#
#     job = dataproc_job(
#         'combine_gvcfs.py',
#         params=dict(
#             cohort_tsv=get_cohort().to_tsv(),
#             out_vds=out_vds_path,
#         ),
#         num_workers=0,
#         autoscaling_policy=f'vcf-combiner-{autoscaling_workers}',
#         long=True,
#     )
#     return job
