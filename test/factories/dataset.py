from typing import Sequence

from cpg_workflows.targets import Cohort, Dataset, SequencingGroup

from .types import DatasetId

DEFAULT_DATASET_NAME: DatasetId = 'local-test'


def create_dataset(
    name: DatasetId = DEFAULT_DATASET_NAME,
    cohort: Cohort | None = None,
    sequencing_groups: Sequence[SequencingGroup] | None = None,
    active: bool = True,
    forced: bool = False,
) -> Dataset:
    """
    Create a dataset with the given name, and add the given sequencing groups to it.
    Also adds the dataset to the Cohort if one is provided.

    Args:
        name (DatasetId, optional):
            The name of your dataset. Make sure it matches the name of the dataset or
            input datasets in your workflow config. Defaults to `'local-test'`.

        cohort (Cohort | None, optional):
            The cohort to add the dataset to. Defaults to `None`.

        sequencing_groups (Sequence[SequencingGroup] | None, optional):
            An iterable of sequencing groups to add to the dataset. Defaults to `None`.

        active (bool, optional):
            Sets to active, so this dataset will be included in the workflow that
            it uses this dataset. Defaults to `True`.

        forced (bool, optional):
            Sets to forced, so this dataset is processed even if outputs exist for it.
            Defaults to `False`.

    Returns:
        Dataset
    """
    dataset = Dataset(name, cohort=cohort)
    if cohort:
        cohort.add_dataset(dataset)

    sequencing_groups = sequencing_groups or []
    for sequencing_group in sequencing_groups:
        dataset._sequencing_group_by_id[sequencing_group.id] = sequencing_group
        sequencing_group.dataset = dataset

    dataset.active = active
    dataset.forced = forced

    return dataset
