"""
Helpers to communicate with the metamist database.
"""

import logging
import pprint
import traceback
from dataclasses import dataclass
from enum import Enum
from typing import Any, Callable, Optional

from gql.transport.exceptions import TransportServerError
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_workflows.filetypes import (
    AlignmentInput,
    BamPath,
    CramPath,
    FastqPair,
    FastqPairs,
)
from cpg_workflows.utils import exists
from metamist import models
from metamist.apis import AnalysisApi
from metamist.exceptions import ApiException, ServiceException
from metamist.graphql import gql, query

GET_SEQUENCING_GROUPS_QUERY = gql(
    """
        query SGQuery($metamist_proj: String!, $only_sgs: [String!]!, $skip_sgs: [String!]!, $sequencing_type: String!) {
            project(name: $metamist_proj) {
                sequencingGroups(id: { in_: $only_sgs, nin: $skip_sgs}, type:  {eq: $sequencing_type}) {
                    id
                    meta
                    platform
                    technology
                    type
                    sample {
                        externalId
                        participant {
                            id
                            externalId
                            phenotypes
                            reportedSex
                            meta
                        }
                    }
                    assays {
                        id
                        meta
                        type
                    }
                }
            }
        }
        """,
)

GET_SEQUENCING_GROUPS_BY_COHORT_QUERY = gql(
    """
    query SGByCohortQuery($cohort_id: String!) {
        cohorts(id: {eq: $cohort_id}) {
            name
            sequencingGroups {
                id
                meta
                platform
                technology
                type
                sample {
                    project {
                        name
                    }
                    externalId
                    participant {
                        id
                        externalId
                        phenotypes
                        reportedSex
                        meta
                    }
                }
                assays {
                    id
                    meta
                    type
                }
            }
        }
    }
    """,
)


GET_ANALYSES_QUERY = gql(
    """
        query AnalysesQuery($metamist_proj: String!, $analysis_type: String!, $analysis_status: AnalysisStatus!) {
            project(name: $metamist_proj) {
                analyses (active: {eq: true}, type: {eq: $analysis_type}, status: {eq: $analysis_status}) {
                    id
                    type
                    meta
                    output
                    status
                    sequencingGroups {
                        id
                    }
                }
            }
        }
        """,
)

GET_PEDIGREE_QUERY = gql(
    """
        query PedigreeQuery($metamist_proj: String!){
            project(name: $metamist_proj) {
                pedigree(replaceWithFamilyExternalIds: false)
            }
        }
    """,
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

    Re-defined in a separate module to decouple from the main metamist module,
    so decorators can use `@stage(analysis_type=AnalysisType.QC)` without importing
    the metamist package.
    """

    QC = 'qc'
    JOINT_CALLING = 'joint-calling'
    GVCF = 'gvcf'
    CRAM = 'cram'
    MITO_CRAM = 'mito-cram'
    CUSTOM = 'custom'
    ES_INDEX = 'es-index'
    COMBINER = 'combiner'

    @staticmethod
    def parse(val: str) -> 'AnalysisType':
        """
        Parse str and create a AnalysisStatus object
        """
        d = {v.value: v for v in AnalysisType}
        if val not in d:
            raise MetamistError(f'Unrecognised analysis type {val}. Available: {list(d.keys())}')
        return d[val.lower()]


@dataclass
class Analysis:
    """
    Metamist DB Analysis entry.

    See the metamist package for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: int
    type: AnalysisType
    status: AnalysisStatus
    sequencing_group_ids: set[str]
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
            sequencing_group_ids=set([s['id'] for s in data['sequencingGroups']]),
            output=output,
            meta=data.get('meta') or {},
        )
        return a


def sort_sgs_by_project(response_data) -> dict:
    """
    Create dictionary organising sequencing groups by project
    {project_id: [sequencing_group_1, sequencing_group_2, ...]}
    """
    result_dict: dict[str, list[str]] = {}

    for sequencing_group in response_data:
        project_id = sequencing_group['sample']['project']['name']

        if project_id not in result_dict:
            result_dict[project_id] = []

        result_dict[project_id].append(sequencing_group)

    return result_dict


class Metamist:
    """
    Communication with metamist.
    """

    def __init__(self) -> None:
        self.default_dataset: str = get_config()['workflow']['dataset']
        self.aapi = AnalysisApi()

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=3, min=8, max=30),
        retry=retry_if_exception_type(ServiceException),
        reraise=True,
    )
    def make_retry_aapi_call(self, api_func: Callable, **kwargv: Any):
        """
        Make a generic API call to self.aapi with retries.
        Retry only if ServiceException is thrown

        TODO: How many retries?
        e.g. try 3 times, wait 2^3: 8, 16, 24 seconds
        """
        try:
            return api_func(**kwargv)
        except ServiceException:
            # raise here so the retry occurs
            logging.warning(
                f'Retrying {api_func} ...',
            )
            raise

    def make_aapi_call(self, api_func: Callable, **kwargv: Any):
        """
        Make a generic API call to self.aapi.
        This is a wrapper around retry of API call to handle exceptions and logging.
        """
        try:
            return self.make_retry_aapi_call(api_func, **kwargv)
        except (ServiceException, ApiException) as e:
            # Metamist API failed even after retries
            # log the error and continue
            traceback.print_exc()
            logging.error(
                f'Error: {e} Call {api_func} failed with payload:\n{str(kwargv)}',
            )
        # TODO: discuss should we catch all here as well?
        # except Exception as e:
        #     # Other exceptions?

        return None

    def get_sgs_for_cohorts(self, cohort_ids: list[str]) -> dict[str, dict[str, Any]]:
        """
        Retrieve the sequencing groups per dataset for a list of cohort IDs.
        """
        return {cohort_id: self.get_sgs_by_project_from_cohort(cohort_id) for cohort_id in cohort_ids}

    def get_sgs_by_project_from_cohort(self, cohort_id: str) -> dict:
        """
        Retrieve sequencing group entries for a cohort.
        """
        entries = query(GET_SEQUENCING_GROUPS_BY_COHORT_QUERY, {'cohort_id': cohort_id})

        # Create dictionary keying sequencing groups by project and including cohort name
        # {
        #     "sequencing_groups": {
        #         project_id: [sequencing_group_1, sequencing_group_2, ...],
        #         ...
        #     },
        #     "name": "CohortName"
        # }
        if len(entries['cohorts']) != 1:
            raise MetamistError('We only support one cohort at a time currently')
        sequencing_groups = entries['cohorts'][0]['sequencingGroups']
        cohort_name = entries['cohorts'][0]['name']
        # TODO (mwelland): future optimisation following closure of #860
        # TODO (mwelland): return all the SequencingGroups in the Cohort, no need for stratification
        return {
            'sequencing_groups': sort_sgs_by_project(sequencing_groups),
            'name': cohort_name,
        }

    @retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=3, min=8, max=30),
        retry=retry_if_exception_type(TransportServerError),
        reraise=True,
    )
    def get_sg_entries(self, dataset_name: str) -> list[dict]:
        """
        Retrieve sequencing group entries for a dataset, in the context of access level
        and filtering options.
        """
        metamist_proj = self.get_metamist_proj(dataset_name)
        logging.info(f'Getting sequencing groups for dataset {metamist_proj}')

        skip_sgs = get_config()['workflow'].get('skip_sgs', [])
        only_sgs = get_config()['workflow'].get('only_sgs', [])
        sequencing_type = get_config()['workflow'].get('sequencing_type')

        if only_sgs and skip_sgs:
            raise MetamistError('Cannot specify both only_sgs and skip_sgs in config')

        sequencing_group_entries = query(
            GET_SEQUENCING_GROUPS_QUERY,
            variables={
                'metamist_proj': metamist_proj,
                'only_sgs': only_sgs,
                'skip_sgs': skip_sgs,
                'sequencing_type': sequencing_type,
            },
        )

        return sequencing_group_entries['project']['sequencingGroups']

    def update_analysis(self, analysis: Analysis, status: AnalysisStatus):
        """
        Update "status" of an Analysis entry.
        """
        self.make_aapi_call(
            self.aapi.update_analysis,
            analysis_id=analysis.id,
            analysis_update_model=models.AnalysisUpdateModel(
                status=models.AnalysisStatus(status.value),
            ),
        )
        # Keeping this as is for compatibility with the existing code
        # However this should only be set after the API call is successful
        analysis.status = status

    # NOTE: This isn't used anywhere.
    def find_joint_calling_analysis(
        self,
        sequencing_group_ids: list[str],
        dataset: str | None = None,
    ) -> Analysis | None:
        """
        Query the DB to find the last completed joint-calling analysis for the sequencing groups.
        """
        metamist_proj = self.get_metamist_proj(dataset)

        data = self.make_aapi_call(
            self.aapi.get_latest_complete_analysis_for_type,
            project=metamist_proj,
            analysis_type=models.AnalysisType('joint-calling'),
        )
        if data is None:
            return None

        a = Analysis.parse(data)
        if not a:
            return None

        assert a.type == AnalysisType.JOINT_CALLING, data
        assert a.status == AnalysisStatus.COMPLETED, data
        if a.sequencing_group_ids != set(sequencing_group_ids):
            return None
        return a

    def get_analyses_by_sgid(
        self,
        sg_ids: list[str],
        analysis_type: AnalysisType,
        analysis_status: AnalysisStatus = AnalysisStatus.COMPLETED,
        meta: dict | None = None,
        dataset: str | None = None,
    ) -> dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type, sequencing group ids,
        and sequencing type, one Analysis object per sequencing group. Assumes the analysis
        is defined for a single sequencing group (that is, analysis_type=cram|gvcf|qc).
        """
        metamist_proj = self.get_metamist_proj(dataset)

        analyses = query(
            GET_ANALYSES_QUERY,
            variables={
                'metamist_proj': metamist_proj,
                'analysis_type': analysis_type.value,
                'analysis_status': analysis_status.name,
            },
        )

        raw_analyses = analyses['project']['analyses']
        if meta:
            raw_analyses = filter_analyses_by_meta(raw_analyses, meta)

        analysis_per_sid: dict[str, Analysis] = {}

        for analysis in raw_analyses:
            a = Analysis.parse(analysis)
            if not a:
                continue

            assert a.status == analysis_status, analysis
            assert a.type == analysis_type, analysis
            if len(a.sequencing_group_ids) < 1:
                logging.warning(f'Analysis has no sequencing group ids. {analysis}')
                continue

            assert len(a.sequencing_group_ids) == 1, analysis
            analysis_per_sid[list(a.sequencing_group_ids)[0]] = a

        logging.info(
            f'Querying {analysis_type} analysis entries for {metamist_proj}: found {len(analysis_per_sid)}',
        )
        return analysis_per_sid

    def create_analysis(
        self,
        output: Path | str,
        type_: str | AnalysisType,
        status: str | AnalysisStatus,
        sequencing_group_ids: list[str],
        dataset: str | None = None,
        meta: dict | None = None,
    ) -> int | None:
        """
        Tries to create an Analysis entry, returns its id if successful.
        """
        metamist_proj = self.get_metamist_proj(dataset)

        if isinstance(type_, AnalysisType):
            type_ = type_.value
        if isinstance(status, AnalysisStatus):
            status = status.value

        am = models.Analysis(
            type=type_,
            status=models.AnalysisStatus(status),
            output=str(output),
            sequencing_group_ids=list(sequencing_group_ids),
            meta=meta or {},
        )
        aid = self.make_aapi_call(
            self.aapi.create_analysis,
            project=metamist_proj,
            analysis=am,
        )
        if aid is None:
            logging.error(
                f'Failed to create Analysis(type={type_}, status={status}, output={str(output)}) in {metamist_proj}',
            )
            return None
        else:
            logging.info(
                f'Created Analysis(id={aid}, type={type_}, status={status}, '
                f'output={str(output)}) in {metamist_proj}',
            )
            return aid

    # NOTE: I don't think this function is used anywhere ~ vivbak 12/06/2023
    def process_existing_analysis(
        self,
        sequencing_group_ids: list[str],
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

        @param sequencing_group_ids: sequencing_group_ids IDs to pull the analysis for
        @param completed_analysis: existing completed analysis of this type for these
        sequencing groups
        @param analysis_type: cram, gvcf, joint_calling
        @param expected_output_fpath: where the workflow expects the analysis output
        file to sit on the bucket (will invalidate the analysis when it doesn't match)
        @param dataset: the name of the dataset to create a new analysis
        @return: path to the output if it can be reused, otherwise None
        """
        label = f'type={analysis_type}'
        if len(sequencing_group_ids) > 1:
            label += f' for {", ".join(sequencing_group_ids)}'

        found_output_fpath: Path | None = None
        if not completed_analysis:
            logging.warning(
                f'Not found completed analysis {label} for '
                f'{f"sequencing group {sequencing_group_ids}" if len(sequencing_group_ids) == 1 else f"{len(sequencing_group_ids)} sequencing_groups" }',
            )
        elif not completed_analysis.output:
            logging.error(f'Found a completed analysis {label}, but the "output" field does not exist or empty')
        else:
            found_output_fpath = completed_analysis.output
            if found_output_fpath != expected_output_fpath:
                logging.error(
                    f'Found a completed analysis {label}, but the "output" path '
                    f'{found_output_fpath} does not match the expected path '
                    f'{expected_output_fpath}',
                )
                found_output_fpath = None
            elif not exists(found_output_fpath):
                logging.error(
                    f'Found a completed analysis {label}, '
                    f'but the "output" file {found_output_fpath} does not exist',
                )
                found_output_fpath = None

        # completed and good exists, can reuse
        if found_output_fpath:
            logging.info(f'Completed analysis {label} exists, reusing the result {found_output_fpath}')
            return found_output_fpath

        # can't reuse, need to invalidate
        if completed_analysis:
            logging.warning(
                f'Invalidating the analysis {label} by setting the status to "failure", '
                f'and resubmitting the analysis.',
            )
            self.update_analysis(completed_analysis, status=AnalysisStatus.FAILED)

        # can reuse, need to create a completed one?
        if exists(expected_output_fpath):
            logging.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed',
            )
            self.create_analysis(
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sequencing_group_ids=sequencing_group_ids,
                dataset=dataset or self.default_dataset,
            )
            return expected_output_fpath

        # proceeding with the standard workflow (creating status=queued, submitting jobs)
        else:
            logging.info(
                f'Expected output file {expected_output_fpath} does not exist, so queueing analysis {label}',
            )
            return None

    def get_ped_entries(self, dataset: str | None = None) -> list[dict[str, str]]:
        """
        Retrieve PED lines for a specified SM project, with external participant IDs.
        """
        metamist_proj = self.get_metamist_proj(dataset)
        entries = query(GET_PEDIGREE_QUERY, variables={'metamist_proj': metamist_proj})

        pedigree_entries = entries['project']['pedigree']

        return pedigree_entries

    def get_metamist_proj(self, dataset: str | None = None) -> str:
        """
        Return the Metamist project name, appending '-test' if the access level is 'test'.
        """
        metamist_proj = dataset or self.default_dataset
        if get_config()['workflow']['access_level'] == 'test' and not metamist_proj.endswith('-test'):
            metamist_proj += '-test'

        return metamist_proj


@dataclass
class Assay:
    """
    Metamist "Assay" entry.

    See metamist for more details:
    https://github.com/populationgenomics/sample-metadata
    """

    id: str
    sequencing_group_id: str
    meta: dict
    assay_type: str
    alignment_input: AlignmentInput | None = None

    @staticmethod
    def parse(
        data: dict,
        sg_id: str,
        check_existence: bool = False,
        run_parse_reads: bool = True,
    ) -> 'Assay':
        """
        Create from a dictionary.
        """

        assay_keys = ['id', 'type', 'meta']
        missing_keys = [key for key in assay_keys if data.get(key) is None]

        if missing_keys:
            raise ValueError(f'Cannot parse metamist Sequence {data}. Missing keys: {missing_keys}')

        assay_type = str(data['type'])
        assert assay_type, data
        mm_seq = Assay(
            id=str(data['id']),
            sequencing_group_id=sg_id,
            meta=data['meta'],
            assay_type=assay_type,
        )
        if run_parse_reads:
            mm_seq.alignment_input = parse_reads(
                sequencing_group_id=sg_id,
                assay_meta=data['meta'],
                check_existence=check_existence,
            )
        return mm_seq


def parse_reads(  # pylint: disable=too-many-return-statements
    sequencing_group_id: str,
    assay_meta: dict,
    check_existence: bool,
) -> AlignmentInput:
    """
    Parse a AlignmentInput object from the meta dictionary.
    `check_existence`: check if fastq/crams exist on buckets.
    Default value is pulled from self.metamist and can be overridden.
    """
    reads_data = assay_meta.get('reads')
    reads_type = assay_meta.get('reads_type')
    reference_assembly = assay_meta.get('reference_assembly', {}).get('location')

    if not reads_data:
        raise MetamistError(f'{sequencing_group_id}: no "meta/reads" field in meta')
    if not reads_type:
        raise MetamistError(f'{sequencing_group_id}: no "meta/reads_type" field in meta')
    supported_types = ('fastq', 'bam', 'cram')
    if reads_type not in supported_types:
        raise MetamistError(
            f'{sequencing_group_id}: ERROR: "reads_type" is expected to be one of {supported_types}',
        )

    if reads_type in ('bam', 'cram'):
        if len(reads_data) > 1:
            raise MetamistError(f'{sequencing_group_id}: supporting only single bam/cram input')

        location = reads_data[0]['location']
        if not (location.endswith('.cram') or location.endswith('.bam')):
            raise MetamistError(
                f'{sequencing_group_id}: ERROR: expected the file to have an extension '
                f'.cram or .bam, got: {location}',
            )
        if get_config()['workflow']['access_level'] == 'test':
            location = location.replace('-main-upload/', '-test-upload/')
        if check_existence and not exists(location):
            raise MetamistError(f'{sequencing_group_id}: ERROR: index file does not exist: {location}')

        # Index:
        index_location = None
        if reads_data[0].get('secondaryFiles'):
            index_location = reads_data[0]['secondaryFiles'][0]['location']
            if (location.endswith('.cram') and not index_location.endswith('.crai')) or (
                location.endswith('.bam') and not index_location.endswith('.bai')
            ):
                raise MetamistError(
                    f'{sequencing_group_id}: ERROR: expected the index file to have an extension '
                    f'.crai or .bai, got: {index_location}',
                )
            if get_config()['workflow']['access_level'] == 'test':
                index_location = index_location.replace('-main-upload/', '-test-upload/')
            if check_existence and not exists(index_location):
                raise MetamistError(f'{sequencing_group_id}: ERROR: index file does not exist: {index_location}')

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
                    f'Sequence data for sequencing group {sequencing_group_id} is incorrectly '
                    f'formatted. Expecting 2 entries per lane (R1 and R2 fastqs), '
                    f'but got {len(lane_pair)}. '
                    f'Read data: {pprint.pformat(lane_pair)}',
                )
            if check_existence and not exists(lane_pair[0]['location']):
                raise MetamistError(
                    f'{sequencing_group_id}: ERROR: read 1 file does not exist: {lane_pair[0]["location"]}',
                )
            if check_existence and not exists(lane_pair[1]['location']):
                raise MetamistError(
                    f'{sequencing_group_id}: ERROR: read 2 file does not exist: {lane_pair[1]["location"]}',
                )

            fastq_pairs.append(
                FastqPair(
                    to_path(lane_pair[0]['location']),
                    to_path(lane_pair[1]['location']),
                ),
            )

        return fastq_pairs


def filter_analyses_by_meta(analyses: list[dict], meta_filter: dict) -> list[dict]:
    """Filter analyses by matching key-value pairs in the meta field."""
    return [a for a in analyses if a.get('meta') and all(a['meta'].get(k) == v for k, v in meta_filter.items())]
