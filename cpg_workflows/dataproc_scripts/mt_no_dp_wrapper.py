from argparse import ArgumentParser

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch


def main():
    parser = ArgumentParser(description='Argument Parser for the ES generation script')
    parser.add_argument('--mt_path', help='MT path name', required=True)
    parser.add_argument('--index', help='ES index name', required=True)
    parser.add_argument('--flag', help='ES index "DONE" file path')
    args = parser.parse_args()

    job = get_batch().new_job(f'Generate {args.index} from {args.mt_path}')
    job.image(config_retrieve(['workflow', 'driver_image']))
    job.storage(config_retrieve(['workflow', 'storage_requirement'], '10Gi'))
    job.memory(config_retrieve(['workflow', 'memory_requirement'], 'highmem'))

    ncpu = config_retrieve(['workflow', 'ncpu'], 4)
    job.cpu(ncpu)

    # localise the MT
    mt_name = args.mt_path.split('/')[-1]
    job.command(f'gcloud --no-user-output-enabled storage cp -r {args.mt_path} $BATCH_TMPDIR')
    job.command(f'mt_to_es --mt_path "${{BATCH_TMPDIR}}/{mt_name}" --index {args.index} --flag {args.flag}')
    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
