import hail as hl

from cpg_workflows.utils import can_reuse


def export_plink(
    dense_mt_path: str,
    output_path: str,
):

    if can_reuse(output_path):
        return []

    dense_mt = hl.read_matrix_table(dense_mt_path)
    return hl.export_plink(dense_mt, output_path, ind_id=dense_mt.s)
