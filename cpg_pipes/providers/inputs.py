"""
Abstract provider of input data sources.
"""

import csv
import logging
from abc import ABC, abstractmethod
from enum import Enum

from ..utils import exists
from ..types import FastqPair, CramPath, AlignmentInput, FastqPairs
from ..targets import Cohort, Dataset, Sex, SequencingType, PedigreeInfo
from .. import Path

logger = logging.getLogger(__file__)


class InputProviderError(Exception):
    """
    Exception thrown when there is something wrong while parsing inputs.
    """

    pass


class InputProvider(ABC):
    """
    Abstract class for implementing input sources.
    """

    def __init__(self, check_files: bool = True):
        self.check_files = check_files

    def populate_cohort(
        self,
        cohort: Cohort,
        dataset_names: list[str] | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        skip_datasets: list[str] | None = None,
        ped_files: list[Path] | None = None,
    ):
        """
        Add datasets in the cohort. There exists only one cohort for
        the pipeline run.
        """
        if dataset_names is not None:
            # Specific datasets requested, so initialising them in advance.
            for ds_name in dataset_names:
                if skip_datasets and ds_name in skip_datasets:
                    logger.info(f'Requested to skipping dataset {ds_name}')
                    continue
                dataset = Dataset.create(ds_name, cohort.namespace)
                entries = self._get_entries(
                    dataset,
                    skip_samples=skip_samples,
                    only_samples=only_samples,
                    skip_datasets=skip_datasets,
                )
                if entries:
                    dataset = cohort.add_dataset(dataset)
                    for entry in entries:
                        self._add_sample(dataset, entry)

        else:
            # We don't know dataset names in advance, so getting all entries.
            entries = self._get_entries(
                skip_samples=skip_samples,
                only_samples=only_samples,
                skip_datasets=skip_datasets,
            )
            for entry in entries:
                ds_name = self.get_dataset_name(entry) or cohort.analysis_dataset.name
                dataset = cohort.create_dataset(ds_name)
                self._add_sample(dataset, entry)

        if not cohort.get_datasets():
            msg = 'No active datasets populated'
            if skip_samples or only_samples or skip_datasets:
                msg += ' (after skipping/picking samples)'
            logger.warning(msg)
            return

        self.populate_alignment_inputs(cohort)
        if cohort.sequencing_type:
            self.filter_sequencing_type(cohort, cohort.sequencing_type)
        self.populate_analysis(cohort)
        self.populate_participants(cohort)
        self.populate_pedigree(cohort)
        if ped_files:
            self.populate_pedigree_from_ped_files(cohort, ped_files)

    @abstractmethod
    def get_entries(
        self,
        dataset: Dataset | None = None,
    ) -> list[dict]:
        """
        Override this method to get a list of data entries (dicts).
        If dataset is not missing, it should be specific to a dataset.
        """

    @abstractmethod
    def get_dataset_name(self, entry: dict) -> str | None:
        """
        Get name of the dataset.
        """

    @abstractmethod
    def get_sample_id(self, entry: dict) -> str:
        """
        Get sample ID.
        """

    @abstractmethod
    def get_external_id(self, entry: dict) -> str | None:
        """
        Get external sample ID.
        """

    @abstractmethod
    def get_participant_id(self, entry: dict) -> str | None:
        """
        Get participant ID. Can be also set later.
        """

    @abstractmethod
    def get_participant_sex(self, entry: dict) -> Sex | None:
        """
        Get participant sex. Can be also set later.
        """

    @abstractmethod
    def get_sample_meta(self, entry: dict) -> dict:
        """
        Get sample metadata. Can be also set later.
        """

    @abstractmethod
    def get_sequencing_type(self, entry: dict) -> SequencingType | None:
        """
        Get sequencing type. Can be also set later.
        """

    @staticmethod
    def filter_sequencing_type(cohort: Cohort, sequencing_type: SequencingType | None):
        """
        Filtering to the samples with only requested sequencing types.
        """
        if not sequencing_type:
            return
        for s in cohort.get_samples():
            if not s.alignment_input_by_seq_type:
                logger.warning(f'{s}: skipping because no sequencing inputs found')
                s.active = False
                continue

            avail_types = list(s.alignment_input_by_seq_type.keys())
            s.alignment_input_by_seq_type = {
                k: v
                for k, v in s.alignment_input_by_seq_type.items()
                if k == sequencing_type
            }
            if not bool(s.alignment_input_by_seq_type):
                logger.warning(
                    f'{s}: skipping because no inputs with data type '
                    f'"{sequencing_type.value}" found in '
                    f'{[k.value for k in avail_types]}'
                )
                s.active = False

    @abstractmethod
    def populate_alignment_inputs(self, cohort: Cohort) -> None:
        """
        Populate alignment inputs for samples.
        """

    @abstractmethod
    def populate_analysis(self, cohort: Cohort) -> None:
        """
        Populate analysis information.
        """

    @abstractmethod
    def populate_participants(self, cohort: Cohort) -> None:
        """
        Populate participant information.
        """

    @abstractmethod
    def populate_pedigree(self, cohort: Cohort) -> None:
        """
        Populate pedigree data (families, sex, relationships).
        """

    @staticmethod
    def populate_pedigree_from_ped_files(
        cohort: Cohort,
        ped_files: list[Path],
    ) -> None:
        """
        Populates pedigree from provided PED files.
        """
        sample_by_participant_id = dict()
        for s in cohort.get_samples():
            sample_by_participant_id[s.participant_id] = s

        for ped_file in ped_files:
            with ped_file.open() as f:
                for line in f:
                    fields = line.strip().split('\t')[:6]
                    fam_id, sam_id, pat_id, mat_id, sex, phenotype = fields
                    if sam_id not in sample_by_participant_id:
                        continue
                    s = sample_by_participant_id[sam_id]
                    s.pedigree = PedigreeInfo(
                        sample=s,
                        fam_id=fam_id,
                        dad=sample_by_participant_id.get(pat_id),
                        mom=sample_by_participant_id.get(mat_id),
                        sex=Sex.parse(sex),
                        phenotype=phenotype or '0',
                    )

    def _get_entries(
        self,
        dataset: Dataset | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        skip_datasets: list[str] | None = None,
    ) -> list[dict]:
        """
        Helper method to get and filter entries.
        """
        entries = self.get_entries(dataset)
        entries = self._filter_samples(
            entries=entries,
            skip_samples=skip_samples,
            only_samples=only_samples,
            skip_datasets=skip_datasets,
        )
        return entries

    def _filter_samples(
        self,
        entries: list[dict[str, str]],
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        skip_datasets: list[str] | None = None,
    ) -> list[dict[str, str]]:
        """
        Apply the only_samples and skip_samples filters.
        """
        filtered_entries = []
        for entry in entries:
            cpgid = self.get_sample_id(entry)
            extid = self.get_external_id(entry)
            ds = self.get_dataset_name(entry)
            if skip_datasets and ds and ds in skip_datasets:
                logger.info(
                    f'Skipping sample (the dataset is skipped): {ds}|{cpgid}|{extid}'
                )
                continue
            if only_samples:
                if cpgid in only_samples or extid in only_samples:
                    logger.info(f'Picking sample: {ds}|{cpgid}|{extid}')
                else:
                    continue
            if skip_samples:
                if cpgid in skip_samples or extid in skip_samples:
                    logger.info(f'Skipping sample: {ds}|{cpgid}|{extid}')
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
    CSV/TSV field names.
    """

    dataset = 'dataset'
    sample = 'sample'
    external_id = 'external_id'
    participant_id = 'participant_id'
    fqs_r1 = 'fqs_r1'
    fqs_r2 = 'fqs_r2'
    cram = 'cram'
    sex = 'sex'
    sequencing_type = 'seq_type'


class CsvInputProvider(InputProvider):
    """
    Input provider that parses data from a CSV/TSV file.
    """

    def __init__(self, fp, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.d = list(csv.DictReader(fp))
        if len(self.d) == 0:
            raise InputProviderError('Empty metadata file')
        if 'sample' not in self.d[0]:
            raise InputProviderError(
                f'Metadata file must have a header and a column "sample", '
                f'got: {self.d[0]}'
            )

    def get_entries(
        self,
        dataset: Dataset | None = None,
    ) -> list[dict]:
        """
        Return list of data entries. Optionally, specific for a dataset.
        """
        entries = self.d
        if dataset:
            entries = [
                e
                for e in entries
                if (
                    not self.get_dataset_name(e)
                    or self.get_dataset_name(e) == dataset.name
                )
            ]
        return entries

    # noinspection PyMethodMayBeStatic
    def get_dataset_name(self, entry: dict) -> str | None:
        """
        Get name of the dataset.
        """
        return entry.get(FieldMap.dataset.value)

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
        reserved_fields = [f.value for f in FieldMap]
        return {k: v for k, v in entry.items() if k not in reserved_fields}

    def get_sequencing_type(self, entry: dict) -> SequencingType:
        """
        Get sequencing type.
        """
        result = SequencingType.parse(entry[FieldMap.sequencing_type.value])
        assert result
        return result

    def populate_analysis(self, cohort: Cohort) -> None:
        """
        Populate Analysis entries.
        """
        pass

    def populate_participants(self, cohort: Cohort) -> None:
        """
        Populate Participant entries.
        """
        pass

    def populate_pedigree(self, cohort: Cohort) -> None:
        """
        Populate pedigree data.
        """
        pass

    def populate_alignment_inputs(self, cohort: Cohort) -> None:
        """
        Populate sequencing inputs for samples.
        """
        data: dict[tuple[str, SequencingType], AlignmentInput] = {}

        for entry in self.get_entries():
            sid = self.get_sample_id(entry)
            seq_type = self.get_sequencing_type(entry)
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
                if len(fqs1) != len(fqs2):
                    raise InputProviderError(
                        'Numbers of fqs_r1 and fqs_r2 values (pipe-separated) '
                        'must match.'
                    )
                if self.check_files:
                    missing_fastqs = [fq for fq in (fqs1 + fqs2) if not exists(fq)]
                    if missing_fastqs:
                        raise InputProviderError(f'FQs {missing_fastqs} does not exist')

                pairs = FastqPairs(
                    [FastqPair(fq1, fq2) for fq1, fq2 in zip(fqs1, fqs2)]
                )
                existing_pairs = data.get((sid, seq_type))
                if existing_pairs:
                    if not isinstance(existing_pairs, FastqPairs):
                        raise InputProviderError(
                            f'Mixed sequencing inputs for sample {id}, '
                            f'type {seq_type.value}: existing: {existing_pairs}, '
                            f'new: {pairs}'
                        )
                    existing_pairs += pairs
                else:
                    data[(sid, seq_type)] = pairs

            elif cram:
                if self.check_files:
                    if not exists(cram):
                        raise InputProviderError(f'CRAM {cram} does not exist')
                if (sid, seq_type) in data:
                    raise InputProviderError(
                        f'Cannot have multiple CRAM/BAM sequencing inputs of the same'
                        f'type for one sample. Sample: {id}, type {seq_type.value}, '
                        f'data: {data[(sid, seq_type)]}, new data: {cram}'
                    )
                data[(sid, seq_type)] = CramPath(cram)

        for sample in cohort.get_samples():
            for (sid, seq_type), alignment_input in data.items():
                if sid == sample.id:
                    sample.alignment_input_by_seq_type[seq_type] = alignment_input
