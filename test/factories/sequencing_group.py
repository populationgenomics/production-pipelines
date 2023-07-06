from pathlib import Path
from typing import Optional


from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs
from cpg_workflows.targets import Dataset, SequencingGroup

from .types import SequencingType


def create_sequencing_group(
    id: str = "CPG01",
    external_id: str = "SAMPLE1",
    participant_id: str | None = None,
    dataset: str | Dataset = "dummy",
    sequencing_type: SequencingType = "genome",
    alignment_input: Optional[FastqPair | FastqPairs | CramPath | BamPath] = None,
) -> SequencingGroup:
    """Creates a new sequencing group with the specified parameters.

    Args:
        id (str):
            The ID of the sequencing group. Defaults to "CPG01".

        external_id (str):
            The external ID of the sequencing group. Defaults to "SAMPLE1".

        participant_id (str, optional):
            The participant ID of the sequencing group. Defaults to None.

        dataset (str or Dataset):
            The dataset associated with the sequencing group. Defaults to "dummy".

        sequencing_type (SequencingType):
            The type of sequencing performed used as the key in the property
            `alignment_input_by_seq_type` on the new instance. Defaults to "genome".

        alignment_input (FastqPair | FastqPairs | CramPath | BamPath, optional):
            The alignment input data. Defaults to None. If None, the property
            `alignment_input_by_seq_type` on the new instance will be None.

    Returns:
        SequencingGroup
    """

    if isinstance(dataset, str):
        dataset = Dataset(dataset)

    return SequencingGroup(
        id=id,
        external_id=external_id,
        participant_id=participant_id,
        dataset=dataset,
        alignment_input_by_seq_type=(
            {sequencing_type: alignment_input} if alignment_input else None
        ),
    )
