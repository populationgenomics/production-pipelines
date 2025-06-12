import gzip
from pathlib import Path

import pytest

from cpg_workflows.jobs.gatk_sv import rename_sv_ids


def write_zippy_bits(tmp_path: Path, lines: list[str], header: list[str] | None = None) -> Path:
    """
    A helper function to write the contents of a VCF file to a string

    Args:
        tmp_path (str): temp root for new temp file
        lines (list[str]): list of lines to write to the file
        header (list[str]): header lines to write to the file, optional
    """

    if header is None:
        header = ['##fileformat=VCFv4.2\n']

    temp_file = tmp_path / 'temp.vcf.gz'
    with gzip.open(temp_file, 'wt') as f:
        f.writelines(header)
        f.writelines(lines)
    return temp_file


def test_rename_sv_ids_1(tmp_path):

    lines = [
        'chr1\t12345\tID\tA\t<DEL>\t.\t.\tEND=56855888\tother\tsections\n',
    ]

    zipfile = write_zippy_bits(tmp_path, lines)
    newpath = tmp_path / 'output.vcf'
    rename_sv_ids(str(zipfile), str(newpath))

    with gzip.open(newpath, 'rt') as f:
        header = next(f)
        assert header == '##fileformat=VCFv4.2\n'

        content = next(f)
        assert content.split('\t')[2] == 'DEL_1-12345-56855888'


def test_rename_sv_ids_2(tmp_path):

    lines = [
        'chr1\t12345\tID\tA\t<DUP>\t.\t.\tEND=56855888\tother\tsections\n',
    ]

    zipfile = write_zippy_bits(tmp_path, lines)
    newpath = tmp_path / 'output.vcf'
    rename_sv_ids(str(zipfile), str(newpath))

    with gzip.open(newpath, 'rt') as f:
        header = next(f)
        assert header == '##fileformat=VCFv4.2\n'

        content = next(f)
        assert content.split('\t')[2] == 'DUP_1-12345-56855888'


def test_rename_sv_ids_3(tmp_path):

    lines = [
        'chr1\t12345\tID\tA\t<BND>\t.\t.\tCHR2=chr5;END2=999\tother\tsections\n',
    ]

    zipfile = write_zippy_bits(tmp_path, lines)
    newpath = tmp_path / 'output.vcf'
    rename_sv_ids(str(zipfile), str(newpath))

    with gzip.open(newpath, 'rt') as f:
        header = next(f)
        assert header == '##fileformat=VCFv4.2\n'

        content = next(f)
        assert content.split('\t')[2] == 'BND_1-12345_5-999'


def test_rename_sv_ids_4(tmp_path):

    lines = [
        'chr1\t12345\tID\tA\t<INS.ALU.ME>\t.\t.\tSVLEN=999\tother\tsections\n',
    ]

    zipfile = write_zippy_bits(tmp_path, lines)
    newpath = tmp_path / 'output.vcf'
    rename_sv_ids(str(zipfile), str(newpath))

    with gzip.open(newpath, 'rt') as f:
        header = next(f)
        assert header == '##fileformat=VCFv4.2\n'

        content = next(f)
        assert content.split('\t')[2] == 'INS.ALU.ME_1-12345_ins999'


def test_rename_sv_ids_5(tmp_path):
    """
    raises an error - no END2 in the INFO field
    and END is missing from the INFO field
    """

    lines = [
        'chr1\t12345\tID\tA\t<BND>\t.\t.\tCHR2=chr5\tother\tsections\n',
    ]

    zipfile = write_zippy_bits(tmp_path, lines)
    newpath = tmp_path / 'output.vcf'
    with pytest.raises(KeyError):
        rename_sv_ids(str(zipfile), str(newpath))
