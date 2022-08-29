import click
import hail as hl


@click.command()
@click.argument('mt_path')
def main(mt_path: str):
    """
    Run VEP on matrix table MT_PATH
    """
    hl.init(default_reference='GRCh38')

    mt = hl.read_matrix_table(mt_path)

    # Filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = mt.filter_rows(mt.alleles[1] != '*')

    mt = hl.vep(mt, config='file:///vep_data/vep-gcloud.json')
    mt.describe()
    mt.show()

    out_path = f'gs://cpg-reference/vep/104.3/dataproc/test/sample.vcf.vep.mt'
    mt.write(out_path)
