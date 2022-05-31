"""
Helpers to communicate with the sample-metadata database.
"""

import logging
import pprint
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
from ...types import FastqPair, CramPath, AlignmentInput, SequencingType, FastqPairs
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

    def find_joint_calling_analysis(
        self,
        sample_ids: list[str],
        project_name: str | None = None,
    ) -> Analysis | None:
        """
        Query the DB to find the last completed joint-calling analysis for the samples.
        """
        try:
            data = self.aapi.get_latest_complete_analysis_for_type(
                project=project_name or self.project_name,
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

    def find_analyses_by_sid(
        self,
        sample_ids: list[str],
        analysis_type: AnalysisType,
        analysis_status: AnalysisStatus = AnalysisStatus.COMPLETED,
        meta: dict | None = None,
        project_name: str | None = None,
    ) -> dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type and samples,
        one Analysis object per sample. Assumes the analysis is defined for a single
        sample (e.g. cram, gvcf).
        """
        project_name = project_name or self.project_name

        analysis_per_sid: dict[str, Analysis] = dict()

        logger.info(
            f'Querying {analysis_type} analysis entries for dataset {project_name}...'
        )
        datas = self.aapi.query_analyses(
            models.AnalysisQueryModel(
                projects=[project_name],
                sample_ids=sample_ids,
                type=models.AnalysisType(analysis_type.value),
                status=models.AnalysisStatus(analysis_status.value),
                meta=meta or {},
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
        logger.info(
            f'Querying {analysis_type} analysis entries for dataset {project_name}: '
            f'found {len(analysis_per_sid)}'
        )
        return analysis_per_sid

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

    def process_existing_analysis(
        self,
        sample_ids: list[str],
        completed_analysis: Analysis | None,
        analysis_type: str,
        expected_output_fpath: Path,
        project_name: str | None = None,
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
        @param expected_output_fpath: where the pipeline expects the analysis output
        file to sit on the bucket (will invalidate the analysis when it doesn't match)
        @param project_name: the name of the project where to create a new analysis
        @return: path to the output if it can be reused, otherwise None
        """
        label = f'type={analysis_type}'
        if len(sample_ids) > 1:
            label += f' for {", ".join(sample_ids)}'

        found_output_fpath: Path | None = None
        if not completed_analysis:
            logger.warning(
                f'Not found completed analysis {label} for '
                f'{f"sample {sample_ids}" if len(sample_ids) == 1 else f"{len(sample_ids)} samples" }'
            )
        elif not completed_analysis.output:
            logger.error(
                f'Found a completed analysis {label}, '
                f'but the "output" field does not exist or empty'
            )
        else:
            found_output_fpath = completed_analysis.output
            if found_output_fpath != expected_output_fpath:
                logger.error(
                    f'Found a completed analysis {label}, but the "output" path '
                    f'{found_output_fpath} does not match the expected path '
                    f'{expected_output_fpath}'
                )
                found_output_fpath = None
            elif not utils.exists(found_output_fpath):
                logger.error(
                    f'Found a completed analysis {label}, '
                    f'but the "output" file {found_output_fpath} does not exist'
                )
                found_output_fpath = None

        # completed and good exists, can reuse
        if found_output_fpath:
            logger.info(
                f'Completed analysis {label} exists, '
                f'reusing the result {found_output_fpath}'
            )
            return found_output_fpath

        # can't reuse, need to invalidate
        if completed_analysis:
            logger.warning(
                f'Invalidating the analysis {label} by setting the status to "failure", '
                f'and resubmitting the analysis.'
            )
            self.update_analysis(completed_analysis, status=AnalysisStatus.FAILED)

        # can reuse, need to create a completed one?
        if utils.exists(expected_output_fpath):
            logger.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed'
            )
            self.create_analysis(
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sample_ids=sample_ids,
                project_name=project_name or self.project_name,
            )
            return expected_output_fpath

        # proceeding with the standard pipeline (creating status=queued, submitting jobs)
        else:
            logger.info(
                f'Expected output file {expected_output_fpath} does not exist, '
                f'so queueing analysis {label}'
            )
            return None

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
            if alignment_input := SmSequence._parse_reads(
                sample_id=sample_id, 
                meta=data['meta'], 
                check_existence=check_existence,
            ):
                sm_seq.alignment_input = alignment_input
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
            fastq_pairs = FastqPairs()
            for lane_pair in reads_data:
                if len(lane_pair) != 2:
                    raise ValueError(
                        f'Sequence data for sample {sample_id} is incorrectly '
                        f'formatted. Expecting 2 entries per lane (R1 and R2 fastqs), '
                        f'but got {len(lane_pair)}. '
                        f'Read data: {pprint.pformat(reads_data)}'
                    )
                if check_existence and not utils.exists(lane_pair[0]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 1 file does not exist: '
                        f'{lane_pair[0]["location"]}'
                    )
                    return None
                if check_existence and not utils.exists(lane_pair[1]['location']):
                    logger.error(
                        f'{sample_id}: ERROR: read 2 file does not exist: '
                        f'{lane_pair[1]["location"]}'
                    )
                    return None

                fastq_pairs.append(
                    FastqPair(
                        to_path(lane_pair[0]['location']),
                        to_path(lane_pair[1]['location']),
                    )
                )

            return fastq_pairs
