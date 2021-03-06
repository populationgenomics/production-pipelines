"""
InputProvider implementation that pulls data from the sample-metadata database.
"""

import logging
import traceback
from collections import defaultdict
from itertools import groupby

from sample_metadata import ApiException

from .smdb import SMDB, SmSequence
from ..inputs import InputProvider, InputProviderError
from ... import Path
from ...targets import Cohort, Sex, PedigreeInfo, Dataset
from ...types import SequencingType

logger = logging.getLogger(__file__)


class SmdbInputProvider(InputProvider):
    """
    InputProvider implementation that pulls data from the sample-metadata database.
    """

    def __init__(
        self,
        db: SMDB,
        check_files: bool = False,
        smdb_errors_are_fatal: bool = True,
    ):
        super().__init__(check_files=check_files)
        self.db = db
        self.smdb_errors_are_fatal = smdb_errors_are_fatal

    def populate_cohort(
        self,
        cohort: Cohort,
        dataset_names: list[str] | None = None,
        skip_samples: list[str] | None = None,
        only_samples: list[str] | None = None,
        skip_datasets: list[str] | None = None,
        ped_files: list[Path] | None = None,
    ) -> Cohort:
        """
        Overriding the superclass method.
        """
        if not dataset_names:
            raise InputProviderError(
                'Dataset must be provided for SmdbInputProvider.populate_cohort()'
            )
        return super().populate_cohort(
            cohort=cohort,
            dataset_names=dataset_names,
            skip_samples=skip_samples,
            only_samples=only_samples,
            skip_datasets=skip_datasets,
            ped_files=ped_files,
        )

    def get_entries(
        self,
        dataset: Dataset | None = None,
    ) -> list[dict]:
        """
        Return list of data entries.
        """
        if dataset is None:
            raise InputProviderError(
                'SmdbInputProvider: dataset must be provided for get_entries()'
            )
        entries = self.db.get_sample_entries(project_name=dataset.name)
        # Adding "dataset" into entries, needed for `self.get_dataset_name()`:
        for e in entries:
            e['dataset'] = dataset.name
        return entries

    def get_dataset_name(self, entry: dict) -> str:
        """
        Get name of the dataset. Not relevant for SMDB because we pull
        specific datasets by their names.
        """
        return entry['dataset']

    def get_sample_id(self, entry: dict) -> str:
        """
        Get sample ID from a sample dict.
        """
        return entry['id'].strip()

    def get_external_id(self, entry: dict) -> str | None:
        """
        Get external sample ID from a sample dict.
        """
        return entry['external_id'].strip()

    def get_participant_id(self, entry: dict) -> str | None:
        """
        Get participant ID from a sample dict.
        """
        return None  # Unknown before all samples are loaded

    def get_participant_sex(self, entry: dict) -> Sex | None:
        """
        Get participant ID from a sample dict.
        """
        return Sex.parse(entry['sex']) if 'sex' in entry else None

    def get_sample_meta(self, entry: dict) -> dict:
        """
        Get sample metadata..
        """
        return entry.get('meta', {})

    def get_sequencing_type(self, entry: dict) -> SequencingType | None:
        """
        Get sequencing type.
        """
        return None  # Unknown before all samples are loaded

    def populate_alignment_inputs(self, cohort: Cohort) -> None:
        """
        Populate sequencing inputs for samples.
        """
        assert cohort.get_sample_ids()
        try:
            found_seqs: list[dict] = self.db.seqapi.get_sequences_by_sample_ids(
                cohort.get_sample_ids()
            )
        except ApiException:
            if self.smdb_errors_are_fatal:
                raise
            else:
                logger.error(
                    'Getting sequencing data from SMDB resulted in an error. '
                    'Continuing without sequencing data because of flag override. '
                    'However, here is the error: '
                )
                traceback.print_exc()
                found_seqs = []
        found_seqs_by_sid = defaultdict(list)
        for found_seq in found_seqs:
            found_seqs_by_sid[found_seq['sample_id']].append(found_seq)
        
        # Checking missing sequences
        if sample_wo_seq := [
            s for s in cohort.get_samples() if s.id not in found_seqs_by_sid
        ]:
            msg = f'No sequencing data found for samples:\n'
            for ds, samples in groupby(sample_wo_seq, key=lambda s: s.dataset.name):
                msg += (
                    f'\t{ds}, {len(list(samples))} samples: '
                    f'{", ".join([s.id for s in samples])}\n'
                )
            if self.smdb_errors_are_fatal:
                raise InputProviderError(msg)
            else:
                msg += '\nContinuing without sequencing data because lenient=true'
                logger.warning(msg)

        for sample in cohort.get_samples():
            for d in found_seqs_by_sid.get(sample.id, []):
                seq = SmSequence.parse(d, self.check_files)
                if seq.alignment_input:
                    if seq.sequencing_type in sample.alignment_input_by_seq_type:
                        raise InputProviderError(
                            f'{sample}: found more than 1 alignment input with '
                            f'sequencing type: {seq.sequencing_type.value}. Check your '
                            f'input provider to make sure there is only one data source '
                            f'of sequencing type per sample.'
                        )
                    sample.alignment_input_by_seq_type[seq.sequencing_type] = \
                        seq.alignment_input

    def populate_analysis(self, cohort: Cohort) -> None:
        """
        Populate Analysis entries.
        """
        pass

    def populate_participants(self, cohort: Cohort) -> None:
        """
        Populate Participant entries.
        """
        for dataset in cohort.get_datasets():
            pid_sid_multi = (
                self.db.papi.get_external_participant_id_to_internal_sample_id(
                    dataset.name
                )
            )
            participant_by_sid = {}
            for group in pid_sid_multi:
                pid = group[0]
                for sid in group[1:]:
                    participant_by_sid[sid] = pid.strip()

            for sample in dataset.get_samples():
                if pid := participant_by_sid.get(sample.id):
                    sample.participant_id = pid

    def populate_pedigree(self, cohort: Cohort) -> None:
        """
        Populate pedigree data for samples.
        """
        sample_by_participant_id = dict()
        for s in cohort.get_samples():
            sample_by_participant_id[s.participant_id] = s

        for dataset in cohort.get_datasets():
            ped_entries = self.db.get_ped_entries(project_name=dataset.name)
            for ped_entry in ped_entries:
                part_id = str(ped_entry['individual_id'])
                if part_id not in sample_by_participant_id.keys():
                    logger.info(
                        f'Participant {part_id} is not found in populated samples '
                        f'and will be skipped'
                    )
                    continue

                s = sample_by_participant_id[part_id]
                maternal_sample = sample_by_participant_id.get(
                    str(ped_entry['maternal_id'])
                )
                paternal_sample = sample_by_participant_id.get(
                    str(ped_entry['paternal_id'])
                )
                s.pedigree = PedigreeInfo(
                    sample=s,
                    fam_id=ped_entry['family_id'],
                    mom=maternal_sample,
                    dad=paternal_sample,
                    sex=Sex.parse(str(ped_entry['sex'])),
                    phenotype=ped_entry['affected'] or '0',
                )

        for dataset in cohort.get_datasets():
            samples_with_ped = [s for s in dataset.get_samples() if s.pedigree]
            logger.info(
                f'{dataset.name}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(dataset.get_samples())}'
            )
