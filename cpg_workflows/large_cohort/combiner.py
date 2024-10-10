"""Runs the combiner.

Inputs:
    File listing sequenging group IDs corresponding to the gvcf files you want to combine
    File listing VDSes in gs:// format, one per line

Raises:
    ValueError: _description_
    ValueError: _description_
    ValueError: _description_
    ValueError: _description_
    ValueError: _description_
    ValueError: _description_

Returns:
    _type_: _description_
"""

from itertools import chain
from typing import TYPE_CHECKING, Any

from graphql import DocumentNode

from hail.vds import VariantDataset, new_combiner, read_vds

from cpg_utils import Path
from cpg_utils.hail_batch import init_batch
from metamist.graphql import gql, query

if TYPE_CHECKING:  # TCH002 https://docs.astral.sh/ruff/rules/typing-only-third-party-import/
    from hail.vds.combiner.variant_dataset_combiner import VariantDatasetCombiner

# Get a cohort. Samples already in a VDS will be filtered out from the gvcf combiner list
# and used as a VDS input
GET_COHORT_QUERY: DocumentNode = gql(
    """
query getCohort($cohort: String!) {
  cohorts(name: {eq: $cohort}) {
    sequencingGroups {
      id
      analyses(type: {eq: "gvcf"}, status: {eq: COMPLETED}, active: {eq: true}) {
        output
      }
    }
  }
}
""",
)

GET_VDS_ANALYSIS_QUERY: DocumentNode = gql(
    """
query getVDSByAnalysisIds($vds_ids: [Int!]!) {
  analyses(id: {in_: $vds_ids}) {
    output
  }
}
""",
)


def _get_samples_from_vds(input_vds: str) -> list[str]:
    return read_vds(input_vds).variant_data.s.collect()


def _parse_metamist_cohort_output(cohort_query_output: dict[str, Any]) -> list[dict[str, Any]]:
    return [sg for sg in cohort_query_output["cohorts"][0]["sequencingGroups"]]


def _parse_metamist_analysis_output(vds_analysis_query_output: dict[str, Any]) -> list[str]:
    return [analysis["output"] for analysis in vds_analysis_query_output["analyses"]]


def run_combiner(
    output_vds_path: Path,
    sequencing_type: str,
    cohort: str,
    tmp_prfx: str,
    existing_vds_ids: int | None = None,
) -> VariantDataset:
    input_vds_samples: list[str] = []
    use_genome_default_intervals: bool = False
    use_exome_default_intervals: bool = False

    existing_vds: list[str] | None = None
    project_gvcf_paths: list[str]

    cohort_query_results: dict[str, Any] = query(
        GET_COHORT_QUERY,
        variables={
            "cohort": cohort,
        },
    )

    init_batch()
    # Only support a single input cohort
    formatted_query_results: list[dict[str, Any]] = _parse_metamist_cohort_output(cohort_query_results)

    if existing_vds_ids:
        vds_query_results: dict[str, Any] = query(
            GET_VDS_ANALYSIS_QUERY,
            variables={
                "vds_ids": existing_vds_ids,
            },
        )
        existing_vds = _parse_metamist_analysis_output(vds_query_results)
        input_vds_samples = list(chain(input_vds_samples, *[_get_samples_from_vds(vds) for vds in existing_vds]))
        # Remove any samples from the query that are already present in a VDS
        formatted_query_results = [sg for sg in formatted_query_results if sg["id"] not in input_vds_samples]

    project_gvcf_paths = [entry["analyses"][0]["output"] for entry in formatted_query_results]

    if sequencing_type == "genome":
        use_genome_default_intervals = True
    elif sequencing_type == "exome":
        use_exome_default_intervals = True

    combiner: VariantDatasetCombiner = new_combiner(
        output_path=str(output_vds_path),
        gvcf_paths=project_gvcf_paths,
        vds_paths=existing_vds,
        reference_genome="GRCh38",
        temp_path=tmp_prfx,
        use_exome_default_intervals=use_exome_default_intervals,
        use_genome_default_intervals=use_genome_default_intervals,
        force=True,
    )

    combiner.run()
    return read_vds(str(output_vds_path))
