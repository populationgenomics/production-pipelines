"""
Generate new intervals from a MatrixTable based on data distribution
A wrapper for the generate_new_intervals.py script
"""

from argparse import ArgumentParser

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch

parser = ArgumentParser()
parser.add_argument('--mt', help='MatrixTable to generate intervals from')
parser.add_argument('--out', help='Output intervals file - path to write the resulting BED in GCP')
parser.add_argument('--meres', help='.interval_list file with intervals to exclude (Centro/Telomeres)')
parser.add_argument('--intervals', help='Number of intervals to generate', type=int)
parser.add_argument('--max_length', help='Max length of an interval', type=int, default=3000000)
args = parser.parse_args()

job = get_batch().new_bash_job(f'Generate approx. {args.intervals} intervals from {args.mt}')
authenticate_cloud_credentials_in_job(job)
job.image(config_retrieve(['workflow', 'driver_image']))

# localise the centromere/telomere file in the batch
meres_in = get_batch().read_input(args.meres)

job.command(
    'new_intervals_from_mt '
    f'--mt {args.mt} '
    f'--out {job.output} '
    f'--meres_file {meres_in} '
    f'--intervals {args.intervals} '
    f'--max_length {args.max_length}',
)
get_batch().write_output(job.output, args.out)
get_batch().run(wait=False)
