import pytest

from cpg_workflows.scripts.long_read_sniffles_vcf_modifier import translate_var_and_sex_to_cn


@pytest.mark.parametrize(
    'contig,var_type,gt,expected',
    (
        ('chr1', 'DUP', '0/0', 2),
        ('chr1', 'DUP', '0/1', 3),
        ('chr1', 'DUP', '1/1', 4),
        ('chr1', 'OTHER', '0/1', 2),
        ('chr1', 'OTHER', '1/1', 2),
        ('chrX', 'DUP', '0/1', 2),
        ('chrX', 'DUP', '0|0', 1),
        ('chrY', 'DUP', '0/1', 2),
        ('chrY', 'DUP', '0/0', 1),
        ('chrM', 'DUP', '0/1', 2),
    ),
)
def test_calculate_cn_male(contig, var_type, gt, expected):
    assert translate_var_and_sex_to_cn(contig, var_type, gt, 1) == expected


@pytest.mark.parametrize(
    'contig,var_type,gt,expected',
    (
        ('chr1', 'DUP', '0/0', 2),
        ('chr1', 'DUP', '0/1', 3),
        ('chr1', 'DUP', '1/1', 4),
        ('chr1', 'OTHER', '0/1', 2),
        ('chr1', 'OTHER', '1/1', 2),
        ('chrX', 'DUP', '0/1', 3),
        ('chrX', 'DUP', '0|0', 2),
        ('chrY', 'DUP', '0/1', 1),
        ('chrY', 'DUP', '0/0', 0),
        ('chrM', 'DUP', '0/1', 2),
    ),
)
def test_calculate_cn_female(contig, var_type, gt, expected):
    assert translate_var_and_sex_to_cn(contig, var_type, gt, 2) == expected
