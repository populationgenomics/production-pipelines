from typing import Literal

from cpg_workflows.filetypes import BamPath, CramPath, FastqPair, FastqPairs, GvcfPath
from cpg_workflows.targets import Dataset, PedigreeInfo, SequencingGroup, Sex

from .dataset import DEFAULT_DATASET_NAME
from .types import (
    DatasetId,
    SequencingGroupExternalId,
    SequencingGroupId,
    SequencingType,
)


def create_sequencing_group(
    id: SequencingGroupId = 'CPGAAAAAA',
    external_id: SequencingGroupExternalId = 'SAMPLE1',
    dataset: DatasetId | Dataset = DEFAULT_DATASET_NAME,
    participant_id: str | None = None,
    meta: dict | None = None,
    sex: Sex | None = None,
    pedigree: PedigreeInfo | None = None,
    # 2024-07-25 mfranklin: I don't think _generally_ the sequencing group should be
    #   limited to the SequencingType[genome, exome] list. Maybe for specific stages,
    #   but not in general for the tests, but I'll leave it here.
    sequencing_type: SequencingType = 'genome',
    sequencing_technology: str = 'short-read',
    sequencing_platform: str = 'illumina',
    alignment_input: FastqPair | FastqPairs | CramPath | BamPath | None = None,
    gvcf: GvcfPath | Literal['auto'] | None = None,
    cram: CramPath | Literal['auto'] | None = None,
    forced: bool = False,
    active: bool = True,
) -> SequencingGroup:
    """
    Creates a new sequencing group with the specified parameters.

    Args:
        id (SequencingGroupId):
            The ID of the sequencing group. Defaults to `CPGAAAAAA`.

        external_id (SequencingGroupExternalId):
            The external ID of the sequencing group. Defaults to `SAMPLE1`.

        dataset (DatasetId | Dataset):
            The dataset associated with the sequencing group. Defaults to `local-test`.

        participant_id (SequencingGroupExternalName | None, optional):
            The participant ID of the sequencing group. Defaults to `None`.

        meta (dict | None, optional):
            Metadata dictionary for the sequencing group. Defaults to `None`.

        sex (Sex | None, optional):
            PED format SEX field for the sequencing group. Defaults to `None`.

        pedigree (PedigreeInfo | None, optional):
            A pedigree information instance describing how this sequencing group
            is related to other sequencing groups in a cohort or dataset. Defaults to
            `None`.

        sequencing_type (SequencingType):
            The sequencing_type of the group of assays. Defaults to 'genome'.

        alignment_input (FastqPair | FastqPairs | CramPath | BamPath | None, optional):
            The alignment input data. Defaults to `None`.

        gvcf (GvcfPath | Literal['auto'] | None, optional):
            Path to a gVCF file path for the sequencing group that will be used to save
            genotyping stage output to. It's meant be used as an alternative to the
            method `make_gvcf_path` on this instance if you want to use a different
            path. If 'auto' is selected, the path is contructed from the parent
            dataset's path and the 'default' storage configuration key. Defaults to
            `None`.

        cram (CramPath | Literal['auto'] | None, optional):
            Path to a CRAM file for the sequencing group that will be used to
            save alignment stage output to. It's meant be used as an alternative to the
            method `make_cram_path` on this instance if you want to use a different
            path. If 'auto' is selected, the path is contructed from the parent
            dataset's path and the 'default' storage configuration key. Defaults to
            `None`.

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
        sequencing_type=sequencing_type,
        sequencing_technology=sequencing_technology,
        sequencing_platform=sequencing_platform,
        dataset=dataset,
        meta=meta,
        sex=sex,
        pedigree=pedigree,
        forced=forced,
        alignment_input=alignment_input,
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
