"""
Helpers to communicate with the sample-metadata database.
"""

import logging
import traceback
from dataclasses import dataclass
from enum import Enum

from sample_metadata import models
from sample_metadata.apis import (
    SampleApi,
    SequenceApi,
    AnalysisApi,
    ParticipantApi,
    FamilyApi,
)
from sample_metadata.exceptions import ApiException

from ... import Path, to_path
from ... import utils
from ...types import FastqPair, CramPath, AlignmentInput, SequencingType
from ..status import AnalysisStatus

logger = logging.getLogger(__file__)


class SmdbError(Exception):
    """
    Raised for problems interacting with sample-metadata database.
    """


class AnalysisType(Enum):
    """
    Corresponds to SMDB Analysis types:
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

    @staticmethod
    def parse(val: str) -> 'AnalysisType':
        """
        Parse str and create a AnalysisStatus object
        """
        d = {v.value: v for v in AnalysisType}
        if val not in d:
            raise SmdbError(
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

    @staticmethod
    def parse(data: dict) -> 'Analysis':
        """
        Parse data to create an Analysis object.
        """
        req_keys = ['id', 'type', 'status']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Analysis" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')

        output = data.get('output')
        if output:
            output = to_path(output)

        a = Analysis(
            id=int(data['id']),
            type=AnalysisType.parse(data['type']),
            status=AnalysisStatus.parse(data['status']),
            sample_ids=set(data.get('sample_ids', [])),
            output=output,
        )
        return a


class SMDB:
    """
    Communication with the SampleMetadata database.
    """

    def __init__(self, project_name: str | None = None):
        """
        @param project_name: default SMDB project name.
        """
        self.sapi = SampleApi()
        self.aapi = AnalysisApi()
        self.seqapi = SequenceApi()
        self.seqapi = SequenceApi()
        self.papi = ParticipantApi()
        self.fapi = FamilyApi()
        self.project_name = project_name

    def get_sample_entries(
        self,
        project_name: str | None = None,
        active: bool = True,
    ) -> list[dict]:
        """
        Get samples in the project as a list of dictionaries.
        """
        project_name = project_name or self.project_name

        logger.debug(f'Finding samples for dataset {project_name}...')
        sample_entries = self.sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                'project_ids': [project_name],
                'active': active,
            }
        )
        logger.info(
            f'Finding samples for project {project_name}: '
            f'found {len(sample_entries)}'
        )
        return sample_entries

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

    def create_analysis(
        self,
        output: Path | str,
        type_: str | AnalysisType,
        status: str | AnalysisStatus,
        sample_ids: list[str],
        project_name: str | None = None,
    ) -> int | None:
        """
        Tries to create an Analysis entry, returns its id if successful.
        """
        project_name = project_name or self.project_name

        if isinstance(type_, AnalysisType):
            type_ = type_.value
        if isinstance(status, AnalysisStatus):
            status = status.value

        am = models.AnalysisModel(
            type=models.AnalysisType(type_),
            status=models.AnalysisStatus(status),
            output=str(output),
            sample_ids=list(sample_ids),
        )
        try:
            aid = self.aapi.create_new_analysis(project=project_name, analysis_model=am)
        except ApiException:
            traceback.print_exc()
            return None
        else:
            logger.info(
                f'Created analysis of type={type_}, status={status} with ID: {aid}'
            )
            return aid

    def get_ped_entries(self, project_name: str | None = None) -> list[dict[str, str]]:
        """
        Retrieve PED lines for a specified SM project, with external participant IDs.
        """
        project_name = project_name or self.project_name

        families = self.fapi.get_families(project_name)
        family_ids = [family['id'] for family in families]
        ped_entries = self.fapi.get_pedigree(
            internal_family_ids=family_ids,
            response_type='json',
            project=project_name,
            replace_with_participant_external_ids=True,
        )
        return ped_entries


@dataclass
class SmSequence:
    """
    Sample-metadata DB "Sequence" entry.

    See sample-metadata for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: str
    sample_id: str
    meta: dict
    sequencing_type: SequencingType
    alignment_input: AlignmentInput | None = None

    @staticmethod
    def parse(data: dict, check_existence: bool) -> 'SmSequence':
        """
        Parse dictionary to create a SmSequence object.
        """
        req_keys = ['id', 'sample_id', 'meta']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Sequence" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')

        sample_id = data['sample_id']
        sm_seq = SmSequence(
            id=data['id'],
            sample_id=sample_id,
            meta=data['meta'],
            sequencing_type=SequencingType.parse(data['type']),
        )
        if data['meta'].get('reads'):
            sm_seq.alignment_input = SmSequence._parse_reads(
                sample_id=sample_id, meta=data['meta'], check_existence=check_existence
            )
        else:
            logger.warning(
                f'{sample_id} sequence: no meta/reads found with FASTQ information'
            )
        return sm_seq

    @staticmethod
    def _parse_reads(  # pylint: disable=too-many-return-statements
        sample_id: str,
        meta: dict,
        check_existence: bool,
    ) -> AlignmentInput | None:
        """
        Parse a AlignmentInput object from the meta dictionary.

        @param check_existence: check if fastq/crams exist on buckets.
        Default value is pulled from self.smdb and can be overridden.
        """
        reads_data = meta.get('reads')
        reads_type = meta.get('reads_type')
        if not reads_data:
            logger.error(f'{sample_id}: no "meta/reads" field in meta')
            return None
        if not reads_type:
            logger.error(f'{sample_id}: no "meta/reads_type" field in meta')
            return None
        supported_types = ('fastq', 'bam', 'cram')
        if reads_type not in supported_types:
            logger.error(
                f'{sample_id}: ERROR: "reads_type" is expected to be one of '
                f'{supported_types}'
            )
            return None

        if reads_type in ('bam', 'cram'):
            if len(reads_data) > 1:
                logger.error(f'{sample_id}: supporting only single bam/cram input')
                return None

            bam_path = reads_data[0]['location']
            if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
                logger.error(
                    f'{sample_id}: ERROR: expected the file to have an extension '
                    f'.cram or .bam, got: {bam_path}'
                )
                return None
            if check_existence and not utils.exists(bam_path):
                logger.error(
                    f'{sample_id}: ERROR: index file does not exist: {bam_path}'
                )
                return None

            # Index:
            if not reads_data[0].get('secondaryFiles'):
                logger.error(
                    f'{sample_id}: ERROR: bam/cram input is expected to have '
                    f'a non-empty list field "secondaryFile" section with indices'
                )
                return None
            index_path = reads_data[0]['secondaryFiles'][0]['location']
            if (
                bam_path.endswith('.cram')
                and not index_path.endswith('.crai')
                or bam_path.endswith('.bai')
                and not index_path.endswith('.bai')
            ):
                logger.error(
                    f'{sample_id}: ERROR: expected the index file to have an extension '
                    f'.crai or .bai, got: {index_path}'
                )
            if check_existence and not utils.exists(index_path):
                logger.error(
                    f'{sample_id}: ERROR: index file does not exist: {index_path}'
                )
                return None

            return CramPath(bam_path, index_path=index_path)

        else:
            fastq_pairs = []
            for lane_data in reads_data:
                assert len(lane_data) == 2, lane_data
                if check_existence and not utils.exists(lane_data[0]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 1 file does not exist: '
                        f'{lane_data[0]["location"]}'
                    )
                    return None
                if check_existence and not utils.exists(lane_data[1]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 2 file does not exist: '
                        f'{lane_data[1]["location"]}'
                    )
                    return None

                fastq_pairs.append(
                    FastqPair(
                        to_path(lane_data[0]['location']),
                        to_path(lane_data[1]['location']),
                    )
                )

            return fastq_pairs
