from typing import Literal, Optional

from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs, GvcfPath
from cpg_workflows.targets import Dataset, PedigreeInfo, SequencingGroup, Sex

from .dataset import DEFAULT_DATASET_NAME
from .types import (
    DatasetId,
    SequencingGroupExternalId,
    SequencingGroupExternalName,
    SequencingGroupId,
    SequencingType,
)


def create_sequencing_group(
    id: SequencingGroupId = 'CPG000001',
    external_id: SequencingGroupExternalId = 'SAMPLE1',
    participant_id: SequencingGroupExternalName | None = None,
    dataset: DatasetId | Dataset = DEFAULT_DATASET_NAME,
    meta: dict | None = None,
    sex: Sex | None = None,
    pedigree: Optional[PedigreeInfo] = None,
    sequencing_type: SequencingType = 'genome',
    alignment_input: Optional[FastqPair | FastqPairs | CramPath | BamPath] = None,
    gvcf: Optional[GvcfPath | Literal['auto']] = None,
    cram: Optional[CramPath | Literal['auto']] = None,
    forced: bool = False,
    active: bool = True,
) -> SequencingGroup:
    """
    Creates a new sequencing group with the specified parameters.

    Args:
        id (SequencingGroupId):
            The ID of the sequencing group. Defaults to `CPG000001`.

        external_id (SequencingGroupExternalId):
            The external ID of the sequencing group. Defaults to `SAMPLE1`.

        participant_id (SequencingGroupExternalName, optional):
            The participant ID of the sequencing group. Defaults to `None`.

        dataset (DatasetId | Dataset):
            The dataset associated with the sequencing group. Defaults to `local-test`.

        meta (dict, optional):
            Metadata dictionary for the sequencing group. Defaults to `None`.

        sex (Sex, optional):
            PED format SEX field for the sequencing group. Defaults to `None`.

        pedigree (PedigreeInfo, optional):
            A pedigree information instance describing how this sequencing group
            is related to other sequencing groups in a cohort or dataset. Defaults to
            `None`.

        sequencing_type (SequencingType):
            The type of sequencing performed used as the key in the property
            `alignment_input_by_seq_type` on the new instance. Defaults to 'genome'.

        alignment_input (FastqPair | FastqPairs | CramPath | BamPath, optional):
            The alignment input data. Defaults to `None`. If `None`, the property
            `alignment_input_by_seq_type` on the new instance will be `None`.

        gvcf (GvcfPath | Literal['auto'], optional):
            Path to a gVCF file path for the sequencing group that will be used to save
            genotyping stage output to. It's meant be used as an alternative to the
            method `make_gvcf_path` on this instance if you want to use a different
            path. If 'auto' is selected, then the path is constructed from using the
            'default' key in storage configuration option for the parent dataset passed
            in via `dataset`. Defaults to `None`.

        cram (CramPath | Literal['auto'], optional):
            Path to a CRAM file for the sequencing group that will be used to
            save alignment stage output to. It's meant be used as an alternative to the
            method `make_cram_path` on this instance if you want to use a different
            path. If 'auto' is selected, then the path is constructed from using the
            'default' key in storage configuration option for the parent dataset passed
            in via `dataset`. Defaults to `None`.

        active (bool, optional):
            Sets to active, so this sequencing group will be included in the workflow
            that it uses this sequencing group. Defaults to `True`.

        forced (bool, optional):
            Sets to forced, so this sequencing group is processed even if outputs exist
            for it. Defaults to `False`.

    Returns:
        SequencingGroup
    """

    if isinstance(dataset, str):
        dataset = Dataset(dataset)

    sg = SequencingGroup(
        id=id,
        external_id=external_id,
        participant_id=participant_id,
        dataset=dataset,
        meta=meta,
        sex=sex,
        pedigree=pedigree,
        forced=forced,
        alignment_input_by_seq_type=(
            {sequencing_type: alignment_input} if alignment_input else None
        ),
    )

    if gvcf == 'auto':
        sg.gvcf = sg.make_gvcf_path()
    else:
        sg.gvcf = gvcf

    if cram == 'auto':
        sg.cram = sg.make_cram_path()
    else:
        sg.cram = cram

    sg.active = active
    sg.forced = forced

    dataset._sequencing_group_by_id[sg.id] = sg

    return sg
