"""
Copy seqr reference data to gs://cpg-reference.
"""
import click
import hail as hl
import logging

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@click.command()
@click.argument('version')
def main(version):
    """On dataproc: copy seqr reference data"""
    hl.init(default_reference='GRCh38')
    ht = hl.read_table(
        'gs://seqr-reference-data/GRCh38/all_reference_data/v2'
        '/combined_reference_data_grch38-2.0.4.ht'
    )
    assert version
    out_bucket = f'gs://cpg-reference/seqr/{version}'
    ht.write(f'{out_bucket}/combined_reference_data_grch38-2.0.4.ht')
    ht = hl.read_table('gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.ht')
    ht.write(f'{out_bucket}/clinvar.GRCh38.ht')


if __name__ == '__main__':
    main()
