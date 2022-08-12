"""
Helpers to communicate with the sample-metadata database.
"""

import logging
import pprint
import traceback
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from cpg_utils.config import get_config

from sample_metadata import models
from sample_metadata.apis import (
    SampleApi,
    SequenceApi,
    AnalysisApi,
    ParticipantApi,
    FamilyApi,
)
from sample_metadata.exceptions import ApiException

from cpg_pipes import Path, to_path
from cpg_pipes import utils
from cpg_pipes.filetypes import FastqPair, CramPath, AlignmentInput, FastqPairs

logger = logging.getLogger(__file__)


class MetamistError(Exception):
    """
    Raised for problems interacting with metamist.
    """


_metamist: Optional['Metamist'] = None


def get_metamist() -> 'Metamist':
    """Return the cohort object, which is a signleton"""
    global _metamist
    if not _metamist:
        _metamist = Metamist()
    return _metamist


class AnalysisStatus(Enum):
    """
    Corresponds to SMDB Analysis statuses:
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
        self.seqapi = SequenceApi()
        self.papi = ParticipantApi()
        self.fapi = FamilyApi()

    def get_sample_entries(
        self,
        dataset: str | None = None,
        active: bool = True,
    ) -> list[dict]:
        """
        Get samples in the project as a list of dictionaries.
        """
        dataset = dataset or self.default_dataset
        logger.debug(f'Finding samples for dataset {dataset}...')
        body = {
            'project_ids': [dataset],
            'active': active,
        }
        sample_entries = self.sapi.get_samples(body_get_samples=body)
        logger.info(
            f'Finding samples for project {dataset}: ' f'found {len(sample_entries)}'
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
        dataset: str | None = None,
    ) -> Analysis | None:
        """
        Query the DB to find the last completed joint-calling analysis for the samples.
        """
        try:
            data = self.aapi.get_latest_complete_analysis_for_type(
                project=dataset or self.default_dataset,
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
        dataset: str | None = None,
    ) -> dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type and samples,
        one Analysis object per sample. Assumes the analysis is defined for a single
        sample (e.g. cram, gvcf).
        """
        dataset = dataset or self.default_dataset

        analysis_per_sid: dict[str, Analysis] = dict()

        logger.info(
            f'Querying {analysis_type} analysis entries for dataset {dataset}...'
        )
        datas = self.aapi.query_analyses(
            models.AnalysisQueryModel(
                projects=[dataset],
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
            f'Querying {analysis_type} analysis entries for dataset {dataset}: '
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
            aid = self.aapi.create_new_analysis(project=dataset, analysis_model=am)
        except ApiException:
            traceback.print_exc()
            return None
        else:
            logger.info(
                f'Created Analysis(id={aid}, type={type_}, status={status}, '
                f'output={str(output)}) in project {dataset}'
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
        @param expected_output_fpath: where the pipeline expects the analysis output
        file to sit on the bucket (will invalidate the analysis when it doesn't match)
        @param dataset: the name of the project where to create a new analysis
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
                dataset=dataset or self.default_dataset,
            )
            return expected_output_fpath

        # proceeding with the standard pipeline (creating status=queued, submitting jobs)
        else:
            logger.info(
                f'Expected output file {expected_output_fpath} does not exist, '
                f'so queueing analysis {label}'
            )
            return None

    def get_ped_entries(self, dataset: str | None = None) -> list[dict[str, str]]:
        """
        Retrieve PED lines for a specified SM project, with external participant IDs.
        """
        dataset = dataset or self.default_dataset

        families = self.fapi.get_families(dataset)
        family_ids = [family['id'] for family in families]
        ped_entries = self.fapi.get_pedigree(
            internal_family_ids=family_ids,
            response_type='json',
            project=dataset,
            replace_with_participant_external_ids=True,
        )

        return ped_entries


@dataclass
class MmSequence:
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
    def parse(data: dict, check_existence: bool = False) -> 'MmSequence':
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
        sequencing_type = data['type']
        assert sequencing_type, data
        sm_seq = MmSequence(
            id=data['id'],
            sample_id=sample_id,
            meta=data['meta'],
            sequencing_type=sequencing_type,
        )
        if data['meta'].get('reads'):
            if alignment_input := MmSequence._parse_reads(
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
        reference_assembly = meta.get('reference_assembly', {}).get('location')

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
            index_path = None
            if reads_data[0].get('secondaryFiles'):
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

            return CramPath(
                bam_path, index_path=index_path, reference_assembly=reference_assembly
            )

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
