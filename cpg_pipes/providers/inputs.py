"""
Abstract provider if input data source.
"""
import csv
import logging
from abc import ABC, abstractmethod
from enum import Enum

from ..utils import exists
from ..types import FastqPair, CramPath, AlignmentInput
from ..targets import Cohort, Dataset, Sex, SequencingType
from .. import Path

logger = logging.getLogger(__file__)


class InputProviderError(Exception):
    """
    Exception thrown when there is something wrong happened parsing
    inputs.
    """

    pass


class InputProvider(ABC):
    """
    Abstract class for implementing inputs source.
    """

    def populate_cohort(
        self,
        cohort: Cohort,
        dataset_names: list[str] | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        ped_files: list[Path] | None = None,
        do_check_seq_existence: bool = False,
    ) -> Cohort:
        """
        Add datasets in the cohort. There exists only one cohort for
        the pipeline run.
        """
        if dataset_names:
            # Specific datasets requested, so initialising them in advance.
            for ds_name in dataset_names:
                dataset = cohort.add_dataset(ds_name)
                entries = self._get_entries(
                    dataset=dataset,
                    skip_samples=skip_samples,
                    only_samples=only_samples,
                )
                for entry in entries:
                    self._add_sample(dataset, entry)

        else:
            # We don't know dataset names in advance, so getting all entries.
            entries = self._get_entries(
                skip_samples=skip_samples, only_samples=only_samples
            )
            for entry in entries:
                ds_name = self.get_dataset_name(cohort, entry)
                dataset = cohort.add_dataset(ds_name)
                self._add_sample(dataset, entry)

        self.populate_alignment_inputs(cohort, do_check_seq_existence)
        self.populate_analysis(cohort)
        self.populate_pedigree(cohort, ped_files)
        return cohort

    @abstractmethod
    def get_entries(self, dataset: Dataset | None = None) -> list[dict]:
        """
        Overide this method to get a list of data entries (dicts).
        If dataset is not missing, it should be specific to a dataset.
        """

    @abstractmethod
    def get_dataset_name(self, cohort: Cohort, entry: dict) -> str:
        """
        Get name of the dataset.
        """

    @abstractmethod
    def get_sample_id(self, entry: dict) -> str:
        """
        Get sample ID from a sample dict.
        """

    @abstractmethod
    def get_external_id(self, entry: dict) -> str | None:
        """
        Get external sample ID from a sample dict.
        """

    @abstractmethod
    def get_participant_id(self, entry: dict) -> str | None:
        """
        Get participant ID from a sample dict.
        """

    @abstractmethod
    def get_participant_sex(self, entry: dict) -> Sex | None:
        """
        Get participant ID from a sample dict.
        """

    @abstractmethod
    def get_sample_meta(self, entry: dict) -> dict:
        """
        Get sample metadata from a sample dict.
        """

    @abstractmethod
    def populate_alignment_inputs(
        self,
        cohort: Cohort,
        do_check_seq_existence: bool = False,
    ):
        """
        Populate sequencing inputs for samples
        """

    def populate_analysis(self, cohort: Cohort) -> None:
        """
        Populate Analysis entries
        """
        pass

    @abstractmethod
    def populate_pedigree(
        self,
        cohort: Cohort,
        ped_files: list[Path] | None = None,
    ) -> None:
        """
        Populate pedigree data from file
        """

    def _get_entries(
        self,
        dataset: Dataset | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
    ) -> list[dict]:
        """
        Helper method to get and filter entries.
        """
        entries = self.get_entries(dataset)
        entries = self._filter_samples(
            entries=entries,
            skip_samples=skip_samples,
            only_samples=only_samples,
        )
        return entries

    def _filter_samples(
        self,
        entries: list[dict[str, str]],
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
    ) -> list[dict[str, str]]:
        """
        Apply the only_samples and skip_samples filters.
        """
        filtered_entries = []
        for entry in entries:
            cpgid = self.get_sample_id(entry)
            extid = self.get_external_id(entry)
            if only_samples:
                if cpgid in only_samples or extid in only_samples:
                    logger.info(f'Picking sample: {cpgid}|{extid}')
                else:
                    continue
            if skip_samples:
                if cpgid in skip_samples or extid in skip_samples:
                    logger.info(f'Skiping sample: {cpgid}|{extid}')
                    continue
            filtered_entries.append(entry)
        return filtered_entries

    def _add_sample(self, dataset: Dataset, entry: dict):
        """
        Helper method to create a sample in a dataset.
        """
        entry = dict(entry)  # making a copy, so methods can call .pop()
        dataset.add_sample(
            id=str(self.get_sample_id(entry)),
            external_id=str(self.get_external_id(entry)),
            participant_id=self.get_participant_id(entry),
            sex=self.get_participant_sex(entry),
            meta=self.get_sample_meta(entry),
        )


class FieldMap(Enum):
    """
    CSV field names.
    """

    dataset = 'dataset'
    sample = 'sample'
    external_id = 'external_id'
    participant_id = 'participant_id'
    fqs_r1 = 'fqs_r1'
    fqs_r2 = 'fqs_r2'
    cram = 'cram'
    sex = 'sex'
    sequencing_type = 'sequencing_type'


class CsvInputProvider(InputProvider):
    """
    Input provider that parses data from a CSV file.
    """

    def __init__(self, fp):
        super().__init__()

        self.d = list(csv.DictReader(fp))
        if len(self.d) == 0:
            raise InputProviderError('Empty metadata TSV')
        if 'sample' not in self.d[0]:
            raise InputProviderError(
                f'Metadata TSV must have a header and a column "sample", '
                f'got: {self.d[0]}'
            )

    def get_entries(self, dataset: Dataset | None = None) -> list[dict]:
        """
        Return list of data entries. Optionally, specific for a dataset.
        """
        entries = self.d
        if dataset:
            entries = [
                e
                for e in entries
                if self.get_dataset_name(dataset.cohort, e) == dataset.name
            ]
        return entries

    # noinspection PyMethodMayBeStatic
    def get_dataset_name(self, cohort: Cohort, entry: dict) -> str:
        """
        Get name of the dataset.
        """
        return entry.get(FieldMap.dataset.value, cohort.analysis_dataset.name)

    # noinspection PyMethodMayBeStatic
    def get_sample_id(self, entry: dict[str, str]) -> str:
        """
        Get sample ID from a sample dict.
        """
        if 'sample' not in entry:
            raise InputProviderError(
                f'Metadata TSV must have a header and a column '
                f'"{FieldMap.sample.value}", got: {entry}'
            )
        return entry[FieldMap.sample.value]

    # noinspection PyMethodMayBeStatic
    def get_external_id(self, entry: dict) -> str | None:
        """
        Get external sample ID from a sample dict.
        """
        return entry.get(FieldMap.external_id.value, None)

    # noinspection PyMethodMayBeStatic
    def get_participant_id(self, entry: dict) -> str | None:
        """
        Get participant ID from a sample dict.
        """
        return entry.get(FieldMap.participant_id.value, None)

    # noinspection PyMethodMayBeStatic
    def get_participant_sex(self, entry: dict) -> Sex:
        """
        Get participant's sex/gender from a sample dict.
        """
        return Sex.parse(entry.get(FieldMap.sex.value, None))

    # noinspection PyMethodMayBeStatic
    def get_sample_meta(self, entry: dict) -> dict:
        """
        Get sample metadata.
        """
        reserverd_fields = [f.value for f in FieldMap]
        return {k: v for k, v in entry.items() if k not in reserverd_fields}

    def populate_analysis(self, cohort: Cohort) -> None:
        """
        Populate Analysis entries
        """
        pass

    def populate_pedigree(
        self,
        cohort: Cohort,
        ped_files: list[Path] | None = None,
    ) -> None:
        """
        Populate pedigree data
        """
        pass

    def populate_alignment_inputs(
        self,
        cohort: Cohort,
        do_check_seq_existence: bool = False,
    ) -> None:
        """
        Populate sequencing inputs for samples.
        """
        d_by_sid: dict[str, AlignmentInput] = {}
        seq_type_by_sid: dict[str, SequencingType] = {}

        for entry in self.get_entries():
            sid = self.get_sample_id(entry)
            fqs1 = [
                f.strip()
                for f in entry.get(FieldMap.fqs_r1.value, '').split('|')
                if f.strip()
            ]
            fqs2 = [
                f.strip()
                for f in entry.get(FieldMap.fqs_r2.value, '').split('|')
                if f.strip()
            ]
            cram = entry.get(FieldMap.cram.value, None)
            if fqs1:
                if len(fqs1) != len(fqs1):
                    raise InputProviderError(
                        'Numbers of fqs_r1 and fqs_r2 values (pipe-separated) '
                        'must match.'
                    )
                if do_check_seq_existence:
                    for fq in fqs1 + fqs2:
                        if not exists(fq):
                            raise InputProviderError(f'FQ {fq} does not exist')
                d_by_sid[sid] = [FastqPair(fq1, fq2) for fq1, fq2 in zip(fqs1, fqs2)]
            elif cram:
                if do_check_seq_existence:
                    if not exists(cram):
                        raise InputProviderError(f'CRAM {cram} does not exist')
                d_by_sid[sid] = CramPath(cram)

            seq_type_by_sid[sid] = SequencingType.parse(
                entry.get(FieldMap.sequencing_type.value, None)
            )

        for sample in cohort.get_samples():
            sample.alignment_input = d_by_sid.get(sample.id)
            if sample.id:
                sample.sequencing_type = seq_type_by_sid[sample.id]
