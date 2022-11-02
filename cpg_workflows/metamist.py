"""
Helpers to communicate with the sample-metadata database.
"""

import logging
import pprint
import traceback
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from sample_metadata import models
from sample_metadata.apis import (
    SampleApi,
    SequenceApi,
    AnalysisApi,
    ParticipantApi,
    FamilyApi,
)
from sample_metadata.exceptions import ApiException

from cpg_utils.config import get_config
from cpg_utils import Path, to_path

from cpg_workflows.utils import exists
from cpg_workflows.filetypes import (
    FastqPair,
    CramPath,
    BamPath,
    AlignmentInput,
    FastqPairs,
)


_metamist: Optional['Metamist'] = None


def get_metamist() -> 'Metamist':
    """Return the cohort object"""
    global _metamist
    if not _metamist:
        _metamist = Metamist()
    return _metamist


class MetamistError(Exception):
    """
    Error while interacting with Metamist.
    """

    pass


class AnalysisStatus(Enum):
    """
    Corresponds to metamist Analysis statuses:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums/analysis.py#L14-L21
    """

    QUEUED = 'queued'
    IN_PROGRESS = 'in-progress'
    FAILED = 'failed'
    COMPLETED = 'completed'
    UNKNOWN = 'unknown'

    @staticmethod
    def parse(name: str) -> 'AnalysisStatus':
        """
        Parse str and create a AnalysisStatus object
        """
        return {v.value: v for v in AnalysisStatus}[name.lower()]


class AnalysisType(Enum):
    """
    Corresponds to metamist Analysis types:
    https://github.com/populationgenomics/sample-metadata/blob/dev/models/enums
    /analysis.py#L4-L11

    Re-defined in a separate module to decouple from the main sample-metadata module,
    so decorators can use `@stage(analysis_type=AnalysisType.QC)` without importing
    the sample-metadata package.
    """

    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    CUSTOM = 'custom'
    ES_INDEX = 'es-index'

    @staticmethod
    def parse(val: str) -> 'AnalysisType':
        """
        Parse str and create a AnalysisStatus object
        """
        d = {v.value: v for v in AnalysisType}
        if val not in d:
            raise MetamistError(
                f'Unrecognised analysis type {val}. Available: {list(d.keys())}'
            )
        return d[val.lower()]


@dataclass
class Analysis:
    """
    Sample metadata DB Analysis entry.

    See the sample-metadata package for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: int
    type: AnalysisType
    status: AnalysisStatus
    sample_ids: set[str]
    output: Path | None
    meta: dict

    @staticmethod
    def parse(data: dict) -> 'Analysis':
        """
        Parse data to create an Analysis object.
        """
        req_keys = ['id', 'type', 'status']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logging.error(f'"Analysis" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse metamist Sequence {data}')

        output = data.get('output')
        if output:
            output = to_path(output)

        a = Analysis(
            id=int(data['id']),
            type=AnalysisType.parse(data['type']),
            status=AnalysisStatus.parse(data['status']),
            sample_ids=set(data.get('sample_ids', [])),
            output=output,
            meta=data.get('meta') or {},
        )
        return a


class Metamist:
    """
    Communication with metamist.
    """

    def __init__(self):
        self.default_dataset: str = get_config()['workflow']['dataset']
        self.sapi = SampleApi()
        self.aapi = AnalysisApi()
        self.seqapi = SequenceApi()
        self.papi = ParticipantApi()
        self.fapi = FamilyApi()

    def get_sample_entries(self, dataset_name: str) -> list[dict]:
        """
        Retrieve sample entries for a dataset, in the context of access level
        and filtering options.
        """
        metamist_proj = dataset_name
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'

        skip_samples = get_config()['workflow'].get('skip_samples', [])
        only_samples = get_config()['workflow'].get('only_samples', [])

        sample_entries = self.sapi.get_samples(
            body_get_samples={'project_ids': [metamist_proj]}
        )
        sample_entries = _filter_sample_entries(
            sample_entries,
            dataset_name,
            skip_samples=skip_samples,
            only_samples=only_samples,
        )
        return sample_entries

    def get_sequence_entries_by_sid(
        self,
        sample_ids: list[str],
        sequencing_type: str,
    ) -> dict[str, list[dict]]:
        """
        Retrieve sample entries for a dataset, in the context of sample IDs
        and sequencing type.
        """
        entries_by_sid = defaultdict(list)

        entries = self.seqapi.get_sequences_by_sample_ids(sample_ids)
        if isinstance(entries, list):
            entries_list: list[dict] = entries
            for entry in entries_list:
                if str(entry['type']) == sequencing_type:
                    entries_by_sid[entry['sample_id']].append(entry)
        else:
            assert isinstance(entries, dict)
            entries_dict: dict[str, list[dict]] = entries
            for sample_id, sample_sequences in entries_dict.items():
                for seq in sample_sequences:
                    if str(seq['type']) == sequencing_type:
                        entries_by_sid[sample_id].append(seq)

        return entries_by_sid

    def get_participant_entries_by_sid(self, dataset_name: str) -> dict[str, dict]:
        """
        Retrieve participant entries for a dataset, in the context of access level.
        """
        metamist_proj = dataset_name
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'

        pid_sid_multi = self.papi.get_external_participant_id_to_internal_sample_id(
            metamist_proj
        )
        sid_by_pid = {}
        for group in pid_sid_multi:
            pid = group[0].strip()
            for sid in group[1:]:
                sid_by_pid[pid] = sid

        entries = self.papi.get_participants(metamist_proj)
        participant_entry_by_sid = {}
        for entry in entries:
            pid = entry['external_id']
            if not (sid := sid_by_pid.get(pid)):
                # This is an expected behaviour: dummy participant entries might be
                # created to fill in the PED data. We should just ignore participants
                # with no associated samples.
                continue
            participant_entry_by_sid[sid] = entry
        return participant_entry_by_sid

    def update_analysis(self, analysis: Analysis, status: AnalysisStatus):
        """
        Update "status" of an Analysis entry.
        """
        try:
            self.aapi.update_analysis_status(
                analysis.id,
                models.AnalysisUpdateModel(status=models.AnalysisStatus(status.value)),
            )
        except ApiException:
            traceback.print_exc()
        analysis.status = status

    def find_joint_calling_analysis(
        self,
        sample_ids: list[str],
        dataset: str | None = None,
    ) -> Analysis | None:
        """
        Query the DB to find the last completed joint-calling analysis for the samples.
        """
        metamist_proj = dataset or self.default_dataset
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'
        try:
            data = self.aapi.get_latest_complete_analysis_for_type(
                project=metamist_proj,
                analysis_type=models.AnalysisType('joint-calling'),
            )
        except ApiException:
            return None
        a = Analysis.parse(data)
        if not a:
            return None
        assert a.type == AnalysisType.JOINT_CALLING, data
        assert a.status == AnalysisStatus.COMPLETED, data
        if a.sample_ids != set(sample_ids):
            return None
        return a

    def get_analyses_by_sid(
        self,
        sample_ids: list[str],
        analysis_type: AnalysisType,
        analysis_status: AnalysisStatus = AnalysisStatus.COMPLETED,
        meta: dict | None = None,
        dataset: str | None = None,
    ) -> dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type, sample ids,
        and sequencing type, one Analysis object per sample. Assumes the analysis
        is defined for a single sample (that is, analysis_type=cram|gvcf|qc).
        """
        dataset = dataset or self.default_dataset
        metamist_proj = dataset or self.default_dataset
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'

        analysis_per_sid: dict[str, Analysis] = dict()

        logging.info(
            f'Querying {analysis_type} analysis entries for {metamist_proj}...'
        )
        meta = meta or {}
        meta['sequencing_type'] = get_config()['workflow']['sequencing_type']

        datas = self.aapi.query_analyses(
            models.AnalysisQueryModel(
                projects=[metamist_proj],
                sample_ids=sample_ids,
                type=models.AnalysisType(analysis_type.value),
                status=models.AnalysisStatus(analysis_status.value),
                meta=meta,
            )
        )

        for data in datas:
            a = Analysis.parse(data)
            if not a:
                continue
            assert a.status == AnalysisStatus.COMPLETED, data
            assert a.type == analysis_type, data
            assert len(a.sample_ids) == 1, data
            analysis_per_sid[list(a.sample_ids)[0]] = a
        logging.info(
            f'Querying {analysis_type} analysis entries for {metamist_proj}: '
            f'found {len(analysis_per_sid)}'
        )
        return analysis_per_sid

    def create_analysis(
        self,
        output: Path | str,
        type_: str | AnalysisType,
        status: str | AnalysisStatus,
        sample_ids: list[str],
        dataset: str | None = None,
        meta: dict | None = None,
    ) -> int | None:
        """
        Tries to create an Analysis entry, returns its id if successful.
        """
        dataset = dataset or self.default_dataset
        metamist_proj = dataset or self.default_dataset
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'

        if isinstance(type_, AnalysisType):
            type_ = type_.value
        if isinstance(status, AnalysisStatus):
            status = status.value

        am = models.AnalysisModel(
            type=models.AnalysisType(type_),
            status=models.AnalysisStatus(status),
            output=str(output),
            sample_ids=list(sample_ids),
            meta=meta or {},
        )
        try:
            aid = self.aapi.create_new_analysis(
                project=metamist_proj, analysis_model=am
            )
        except ApiException:
            traceback.print_exc()
            return None
        else:
            logging.info(
                f'Created Analysis(id={aid}, type={type_}, status={status}, '
                f'output={str(output)}) in {metamist_proj}'
            )
            return aid

    def process_existing_analysis(
        self,
        sample_ids: list[str],
        completed_analysis: Analysis | None,
        analysis_type: str,
        expected_output_fpath: Path,
        dataset: str | None = None,
    ) -> Path | None:
        """
        Checks whether existing analysis exists, and output matches the expected output
        file. Invalidates bad analysis by setting status=failure, and submits a
        status=completed analysis if the expected output already exists.

        Returns the path to the output if it can be reused, otherwise None.

        @param sample_ids: sample IDs to pull the analysis for
        @param completed_analysis: existing completed analysis of this type for these
        samples
        @param analysis_type: cram, gvcf, joint_calling
        @param expected_output_fpath: where the workflow expects the analysis output
        file to sit on the bucket (will invalidate the analysis when it doesn't match)
        @param dataset: the name of the dataset to create a new analysis
        @return: path to the output if it can be reused, otherwise None
        """
        label = f'type={analysis_type}'
        if len(sample_ids) > 1:
            label += f' for {", ".join(sample_ids)}'

        found_output_fpath: Path | None = None
        if not completed_analysis:
            logging.warning(
                f'Not found completed analysis {label} for '
                f'{f"sample {sample_ids}" if len(sample_ids) == 1 else f"{len(sample_ids)} samples" }'
            )
        elif not completed_analysis.output:
            logging.error(
                f'Found a completed analysis {label}, '
                f'but the "output" field does not exist or empty'
            )
        else:
            found_output_fpath = completed_analysis.output
            if found_output_fpath != expected_output_fpath:
                logging.error(
                    f'Found a completed analysis {label}, but the "output" path '
                    f'{found_output_fpath} does not match the expected path '
                    f'{expected_output_fpath}'
                )
                found_output_fpath = None
            elif not exists(found_output_fpath):
                logging.error(
                    f'Found a completed analysis {label}, '
                    f'but the "output" file {found_output_fpath} does not exist'
                )
                found_output_fpath = None

        # completed and good exists, can reuse
        if found_output_fpath:
            logging.info(
                f'Completed analysis {label} exists, '
                f'reusing the result {found_output_fpath}'
            )
            return found_output_fpath

        # can't reuse, need to invalidate
        if completed_analysis:
            logging.warning(
                f'Invalidating the analysis {label} by setting the status to "failure", '
                f'and resubmitting the analysis.'
            )
            self.update_analysis(completed_analysis, status=AnalysisStatus.FAILED)

        # can reuse, need to create a completed one?
        if exists(expected_output_fpath):
            logging.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed'
            )
            self.create_analysis(
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sample_ids=sample_ids,
                dataset=dataset or self.default_dataset,
            )
            return expected_output_fpath

        # proceeding with the standard workflow (creating status=queued, submitting jobs)
        else:
            logging.info(
                f'Expected output file {expected_output_fpath} does not exist, '
                f'so queueing analysis {label}'
            )
            return None

    def get_ped_entries(self, dataset: str | None = None) -> list[dict[str, str]]:
        """
        Retrieve PED lines for a specified SM project, with external participant IDs.
        """
        metamist_proj = dataset or self.default_dataset
        if get_config()['workflow']['access_level'] == 'test':
            metamist_proj += '-test'

        entries = self.fapi.get_families(metamist_proj)
        family_ids = [entry['id'] for entry in entries]

        # Since `fapi.get_pedigree` is a GET endpoint, it is limited by the length of
        # the request string. It would stall with the number of families above ~600.
        # To mitigate this, we split the input into chunks. 500 families should be
        # a safe number of families in one chunk.
        def _chunks(seq, size):
            return (seq[pos : pos + size] for pos in range(0, len(seq), size))

        ped_entries = []
        chunk_size = 500
        for i, fam_ids_chunk in enumerate(_chunks(family_ids, chunk_size)):
            logging.info(
                f'Running fapi.get_pedigree on families #{i * chunk_size + 1}..'
                f'{i * chunk_size + 1 + len(fam_ids_chunk) - 1} '
                f'(out of {len(family_ids)})'
            )
            ped_entries.extend(
                self.fapi.get_pedigree(
                    internal_family_ids=fam_ids_chunk,
                    export_type='json',
                    project=metamist_proj,
                    replace_with_participant_external_ids=True,
                )
            )

        return ped_entries


@dataclass
class Sequence:
    """
    Metamist "Sequence" entry.

    See metamist for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: str
    sample_id: str
    meta: dict
    sequencing_type: str
    alignment_input: AlignmentInput | None = None

    @staticmethod
    def parse(
        data: dict,
        check_existence: bool = False,
        parse_reads: bool = True,
    ) -> 'Sequence':
        """
        Create from a dictionary.
        """
        req_keys = ['id', 'sample_id', 'meta']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logging.error(f'"Sequence" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse metamist Sequence {data}')

        sample_id = str(data['sample_id'])
        sequencing_type = str(data['type'])
        assert sequencing_type, data
        mm_seq = Sequence(
            id=str(data['id']),
            sample_id=sample_id,
            meta=data['meta'],
            sequencing_type=sequencing_type,
        )
        if parse_reads:
            mm_seq.alignment_input = Sequence.parse_reads(
                sample_id=sample_id,
                meta=data['meta'],
                check_existence=check_existence,
            )
        return mm_seq

    @staticmethod
    def parse_reads(  # pylint: disable=too-many-return-statements
        sample_id: str,
        meta: dict,
        check_existence: bool,
    ) -> AlignmentInput:
        """
        Parse a AlignmentInput object from the meta dictionary.
        `check_existence`: check if fastq/crams exist on buckets.
        Default value is pulled from self.metamist and can be overridden.
        """
        reads_data = meta.get('reads')
        reads_type = meta.get('reads_type')
        reference_assembly = meta.get('reference_assembly', {}).get('location')

        if not reads_data:
            raise MetamistError(f'{sample_id}: no "meta/reads" field in meta')
        if not reads_type:
            raise MetamistError(f'{sample_id}: no "meta/reads_type" field in meta')
        supported_types = ('fastq', 'bam', 'cram')
        if reads_type not in supported_types:
            raise MetamistError(
                f'{sample_id}: ERROR: "reads_type" is expected to be one of '
                f'{supported_types}'
            )

        if reads_type in ('bam', 'cram'):
            if len(reads_data) > 1:
                raise MetamistError(
                    f'{sample_id}: supporting only single bam/cram input'
                )

            location = reads_data[0]['location']
            if not (location.endswith('.cram') or location.endswith('.bam')):
                raise MetamistError(
                    f'{sample_id}: ERROR: expected the file to have an extension '
                    f'.cram or .bam, got: {location}'
                )
            if get_config()['workflow']['access_level'] == 'test':
                location = location.replace('-main-upload/', '-test-upload/')
            if check_existence and not exists(location):
                raise MetamistError(
                    f'{sample_id}: ERROR: index file does not exist: {location}'
                )

            # Index:
            index_location = None
            if reads_data[0].get('secondaryFiles'):
                index_location = reads_data[0]['secondaryFiles'][0]['location']
                if (
                    location.endswith('.cram')
                    and not index_location.endswith('.crai')
                    or location.endswith('.bai')
                    and not index_location.endswith('.bai')
                ):
                    raise MetamistError(
                        f'{sample_id}: ERROR: expected the index file to have an extension '
                        f'.crai or .bai, got: {index_location}'
                    )
                if get_config()['workflow']['access_level'] == 'test':
                    index_location = index_location.replace(
                        '-main-upload/', '-test-upload/'
                    )
                if check_existence and not exists(index_location):
                    raise MetamistError(
                        f'{sample_id}: ERROR: index file does not exist: {index_location}'
                    )

            if location.endswith('.cram'):
                return CramPath(
                    location,
                    index_path=index_location,
                    reference_assembly=reference_assembly,
                )
            else:
                assert location.endswith('.bam')
                return BamPath(location, index_path=index_location)

        else:
            fastq_pairs = FastqPairs()
            for lane_pair in reads_data:
                if len(lane_pair) != 2:
                    raise ValueError(
                        f'Sequence data for sample {sample_id} is incorrectly '
                        f'formatted. Expecting 2 entries per lane (R1 and R2 fastqs), '
                        f'but got {len(lane_pair)}. '
                        f'Read data: {pprint.pformat(reads_data)}'
                    )
                if get_config()['workflow']['access_level'] == 'test':
                    lane_pair[0]['location'] = lane_pair[0]['location'].replace(
                        '-main-upload/', '-test-upload/'
                    )
                    lane_pair[1]['location'] = lane_pair[1]['location'].replace(
                        '-main-upload/', '-test-upload/'
                    )
                if check_existence and not exists(lane_pair[0]['location']):
                    raise MetamistError(
                        f'{sample_id}: ERROR: read 1 file does not exist: '
                        f'{lane_pair[0]["location"]}'
                    )
                if check_existence and not exists(lane_pair[1]['location']):
                    raise MetamistError(
                        f'{sample_id}: ERROR: read 2 file does not exist: '
                        f'{lane_pair[1]["location"]}'
                    )

                fastq_pairs.append(
                    FastqPair(
                        to_path(lane_pair[0]['location']),
                        to_path(lane_pair[1]['location']),
                    )
                )

            return fastq_pairs


def _filter_sample_entries(
    entries: list[dict[str, str]],
    dataset_name: str,
    skip_samples: list[str] | None = None,
    only_samples: list[str] | None = None,
) -> list[dict]:
    """
    Apply the only_samples and skip_samples filters to sample entries.
    """

    filtered_entries = []
    for entry in entries:
        cpgid = entry['id']
        extid = entry['external_id']
        if only_samples:
            if cpgid in only_samples or extid in only_samples:
                logging.info(f'Picking sample: {dataset_name}|{cpgid}|{extid}')
            else:
                continue
        if skip_samples:
            if cpgid in skip_samples or extid in skip_samples:
                logging.info(f'Skipping sample: {dataset_name}|{cpgid}|{extid}')
                continue
        filtered_entries.append(entry)
    return filtered_entries
