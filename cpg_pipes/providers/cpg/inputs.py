"""
InputProvider implementation that pulls data from the sample-metadata database.
"""

import logging
import traceback
from functools import lru_cache

from sample_metadata import ApiException

from .smdb import SMDB, SmSequence
from ..inputs import InputProvider, InputProviderError
from ...targets import Cohort, Sex, PedigreeInfo

logger = logging.getLogger(__file__)


class SmdbInputProvider(InputProvider):
    """
    InputProvider implementation that pulls data from the sample-metadata database.
    """

    def __init__(self, db: SMDB):
        super().__init__()
        self.db = db

    def get_entries(self, dataset_name: str | None = None) -> list[dict]:
        """
        Return list of data entries.
        """
        if dataset_name is None:
            raise InputProviderError(
                'SmdbInputProvider: dataset_name must be provided for get_entries()'
            )
        entries = self.db.get_sample_entries(dataset_name=dataset_name)
        # Adding "dataset" into entries, needed for `self.get_dataset_name()`:
        for e in entries:
            e['dataset'] = dataset_name
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
    
    @lru_cache
    def _get_participant_id_map(self, dataset_name: str) -> dict[str, str]:
        """
        Returns map of participant IDs to internal CPG IDs.
        """
        pid_sid_multi = self.db.papi.get_external_participant_id_to_internal_sample_id(
            dataset_name
        )
        sid_to_pid = {}
        for group in pid_sid_multi:
            pid = group[0]
            for sid in group[1:]:
                sid_to_pid[sid] = pid
        return sid_to_pid

    def get_participant_id(self, entry: dict) -> str | None:
        """
        Get participant ID from a sample dict.
        """
        sid_to_pid = self._get_participant_id_map(self.get_dataset_name(entry))
        pid = sid_to_pid.get(self.get_sample_id(entry))
        return pid.strip() if pid else None

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

    def populate_alignment_inputs(
        self,
        cohort: Cohort,
        do_check_seq_existence: bool = False,
    ) -> None:
        """
        Populate sequencing inputs for samples.
        """
        try:
            seq_infos: list[dict] = self.db.seqapi.get_sequences_by_sample_ids(
                [s.id for s in cohort.get_samples()]
            )
        except ApiException:
            traceback.print_exc()
            return

        seq_info_by_sid = {seq['sample_id']: seq for seq in seq_infos}
        for sample in cohort.get_samples():
            seq_info = seq_info_by_sid[sample.id]
            seq = SmSequence.parse(seq_info, do_check_seq_existence)
            sample.alignment_input = seq.alignment_input
            sample.sequencing_type = seq.sequencing_type

    def populate_pedigree(
        self,
        cohort: Cohort,
    ) -> None:
        """
        Populate pedigree data for samples.
        """
        sample_by_participant_id = dict()
        for s in cohort.get_samples():
            sample_by_participant_id[s.participant_id] = s

        sample_by_internal_id = dict()
        for s in cohort.get_samples():
            sample_by_internal_id[s.id] = s

        for dataset in cohort.get_datasets():
            ped_entries = self.db.get_ped_entries(dataset_name=dataset.name)
            for entry in ped_entries:
                sam_id = entry['individual_id']
                if sam_id not in sample_by_internal_id:
                    continue
                s = sample_by_internal_id[sam_id]
                s.pedigree = PedigreeInfo(
                    sample=s,
                    fam_id=entry['family_id'],
                    mom=sample_by_internal_id.get(entry['maternal_id']),
                    dad=sample_by_internal_id.get(entry['paternal_id']),
                    sex=Sex.parse(str(entry['sex'])),
                    phenotype=entry['affected'] or '0',
                )

        for dataset in cohort.get_datasets():
            samples_with_ped = [s for s in dataset.get_samples() if s.pedigree]
            logger.info(
                f'{dataset.name}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(dataset.get_samples())}'
            )
