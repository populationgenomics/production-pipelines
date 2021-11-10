"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
from dataclasses import dataclass
from typing import List, Dict, Optional, Set, Collection

from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_production_pipelines import utils, resources
from cpg_production_pipelines.jobs import wrap_command
from cpg_production_pipelines.hailbatch import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)



@dataclass
class Analysis:
    """
    Represents the Analysis SampleMetadata DB entry
    """

    id: str
    type: str
    status: str
    sample_ids: Set[str]
    output: Optional[str]

    def make_analysis_model(self):
        from sample_metadata import AnalysisModel
        return AnalysisModel(
            type=self.type,
            output=self.output,
            status=self.status,
            sample_ids=self.sample_ids,
        )


class Sequence:
    """
    Represents the Sequence SampleMetadata DB entry
    """
    
    def __init__(self, id, meta, smdb):
        self.id = id
        self.meta = meta
        self.smdb = smdb
    
    @staticmethod
    def parse(data: Dict, smdb):
        return Sequence(data['id'], data['meta'], smdb)

    def parse_reads_from_metadata(self):
        return parse_reads_from_metadata(
            self.meta,
            check_existence=self.smdb.do_check_existence
        )


class SMDB:
    """
    Abstracting the communication with the SampleMetadata database.
    """

    def __init__(
        self, 
        analysis_project: str, 
        do_update_analyses: bool = True,
        do_check_existence: bool = True,
    ):
        """
        :param analysis_project: project where to create the "analysis" entries.
        :param do_update_analyses: if not set, won't update "analysis" entries' 
            statuses.
        :param do_check_existence: when querying "sequence" or "analysis" entries
            with files, check those files existence with gsutil. For "sequence", will
            throw an error. For "analysis", will invalidate by setting status=failure.
        """
        from sample_metadata import AnalysisApi, SequenceApi, SampleApi
        self.sapi = SampleApi()
        self.aapi = AnalysisApi()
        self.seqapi = SequenceApi()
        self.analysis_project = analysis_project
        self.do_update_analyses = do_update_analyses
        self.do_check_existence = do_check_existence

    def get_samples_by_project(
        self,
        project_names: List[str],
        namespace: str,
        skip_samples: Optional[List[str]] = None,
    ) -> Dict[str, List[Dict]]:
        """
        Returns a dictionary of samples per input projects
        """
        samples_by_project: Dict[str, List[Dict]] = dict()
        for proj_name in project_names:
            logger.info(f'Finding samples for project {proj_name}')
            input_proj_name = proj_name
            if namespace != 'main':
                input_proj_name += '-test'
            samples = self.sapi.get_samples(
                body_get_samples_by_criteria_api_v1_sample_post={
                    'project_ids': [input_proj_name],
                    'active': True,
                }
            )
            samples_by_project[proj_name] = []
            for s in samples:
                if skip_samples and s['id'] in skip_samples:
                    logger.info(f'Skiping sample: {s["id"]}')
                    continue
                samples_by_project[proj_name].append(s)
        return samples_by_project

    def find_seq_info_by_sid(self, sample_ids) -> Dict[List, Sequence]:
        """
        Return a dict of "Sequence" entries by sample ID
        """
        seq_infos: List[Dict] = self.seqapi.get_sequences_by_sample_ids(sample_ids)
        seq_infos = [Sequence.parse(d, self) for d in seq_infos]
        seq_info_by_sid = dict(zip(sample_ids, seq_infos))
        return seq_info_by_sid

    def update_analysis(self, analysis: Analysis, status: str):
        """
        Update "status" of an Analysis entry
        """
        from sample_metadata import AnalysisUpdateModel

        if not self.do_update_analyses:
            return
        self.aapi.update_analysis_status(
            analysis.id, AnalysisUpdateModel(status=status)
        )
        analysis.status = status

    def find_joint_calling_analysis(
        self,
        sample_ids: Collection[str],
    ) -> Optional[Analysis]:
        """
        Query the DB to find the last completed joint-calling analysis for the samples
        """
        data = self.aapi.get_latest_complete_analysis_for_type(
            project=self.analysis_project,
            analysis_type='joint-calling',
        )
        a = _parse_analysis(data)
        if not a:
            return None
        assert a.type == 'joint-calling', data
        assert a.status == 'completed', data
        if a.sample_ids != set(sample_ids):
            return None
        return a

    def find_analyses_by_sid(
        self,
        sample_ids: Collection[str],
        analysis_type: str,
        project: Optional[str] = None,
    ) -> Dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type and samples,
        one Analysis object per sample. Assumes the analysis is defined for a single
        sample (e.g. cram, gvcf)
        """
        from sample_metadata import exceptions
        
        project = project or self.analysis_project
        analysis_per_sid: Dict[str, Analysis] = dict()
        try:
            logger.info(
                f'Querying analysis entries for project {project}'
            )
            datas = self.aapi.get_latest_analysis_for_samples_and_type(
                project=project,
                analysis_type=analysis_type,
                request_body=sample_ids,
            )
        except exceptions.ApiException:
            return dict()

        for data in datas:
            a = _parse_analysis(data)
            if not a:
                continue
            if a.status == 'completed':
                assert a.type == analysis_type, data
                assert len(a.sample_ids) == 1, data
                analysis_per_sid[list(a.sample_ids)[0]] = a
        return analysis_per_sid

    def make_sm_in_progress_job(self, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "in-progress"
        """
        kwargs['status'] = 'in-progress'
        return self.make_sm_update_status_job(*args, **kwargs)

    def make_sm_completed_job(self, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "completed"
        """
        kwargs['status'] = 'completed'
        return self.make_sm_update_status_job(*args, **kwargs)

    def make_sm_update_status_job(
        self,
        b: Batch,
        analysis_id: str,
        analysis_type: str,
        status: str,
        sample_name: Optional[str] = None,
        project_name: Optional[str] = None,
    ) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status.
        """
        assert status in ['in-progress', 'failed', 'completed', 'queued']
        job_name = ''
        if project_name and sample_name:
            job_name += f'{project_name}/{sample_name}: '
        job_name += f'Update SM: {analysis_type} to {status}'

        if not self.do_update_analyses:
            return b.new_job(f'{job_name} [skip]')

        j = b.new_job(job_name)
        j.image(resources.SM_IMAGE)
        j.command(wrap_command(f"""\
        export SM_DEV_DB_PROJECT={self.analysis_project}
        export SM_ENVIRONMENT=PRODUCTION
        
        cat <<EOT >> update.py
        from sample_metadata.api import AnalysisApi
        from sample_metadata import AnalysisUpdateModel
        from sample_metadata import exceptions
        import traceback
        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id='{analysis_id}',
                analysis_update_model=AnalysisUpdateModel(status='{status}'),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        EOT
        python update.py
        """, setup_gcp=True))
        return j

    def create_analysis(
        self,
        type_: str,
        output: str,
        status: str,
        sample_ids: Collection[str],
    ) -> Optional[int]:
        """
        Tries to create an Analysis entry, returns its id if successfuly
        """
        from sample_metadata import AnalysisModel

        if not self.do_update_analyses:
            return None

        am = AnalysisModel(
            type=type_,
            output=output,
            status=status,
            sample_ids=sample_ids,
        )
        aid = self.aapi.create_new_analysis(
            project=self.analysis_project, analysis_model=am
        )
        logger.info(f'Queueing joint-calling with analysis ID: {aid}')
        return aid

    def process_existing_analysis(
        self,
        sample_ids: Collection[str],
        completed_analysis: Optional[Analysis],
        analysis_type: str,
        expected_output_fpath: str,
    ) -> Optional[str]:
        """
        Checks whether existing analysis exists, and output matches the expected output
        file. Invalidates bad analysis by setting status=failure, and submits a
        status=completed analysis if the expected output already exists.

        Returns the path to the output if it can be reused, otherwise None.

        :param sample_ids: sample IDs to pull the analysis for
        :param completed_analysis: existing completed analysis of this type for these samples
        :param analysis_type: cram, gvcf, joint_calling
        :param expected_output_fpath: where the pipeline expects the analysis output file
            to sit on the bucket (will invalidate the analysis if it doesn't match)
        :return: path to the output if it can be reused, otherwise None
        """
        label = f'type={analysis_type}'
        if len(sample_ids) > 1:
            label += f' for {", ".join(sample_ids)}'

        found_output_fpath = None
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
            found_output_fpath = str(completed_analysis.output)
            if found_output_fpath != expected_output_fpath:
                logger.error(
                    f'Found a completed analysis {label}, but the "output" path '
                    f'{found_output_fpath} does not match the expected path '
                    f'{expected_output_fpath}'
                )
                found_output_fpath = None
            elif not utils.file_exists(found_output_fpath):
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
            self.update_analysis(completed_analysis, status='failed')

        # can reuse, need to create a completed one?
        if utils.file_exists(expected_output_fpath):
            logger.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed'
            )
            self.create_analysis(
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sample_ids=sample_ids,
            )
            return expected_output_fpath

        # proceeding with the standard pipeline (creating status=queued, submitting jobs)
        else:
            logger.info(
                f'Expected output file {expected_output_fpath} does not exist, '
                f'so queueing analysis {label}'
            )
            return None

    def add_running_and_completed_update_jobs(
        self,
        b,
        analysis_type,
        output_path,
        sample_names,
        project_name,
        first_j,
        last_j,
        depends_on,
    ):
        if not self.do_update_analyses:
            return last_j
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = self.create_analysis(
            type_=analysis_type,
            output=output_path,
            status='queued',
            sample_ids=sample_names,
        )
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = self.make_sm_in_progress_job(
            b,
            analysis_id=aid,
            analysis_type=analysis_type,
            project_name=project_name,
            sample_name=sample_names[0] if len(sample_names) == 1 else None,
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = self.make_sm_completed_job(
            b,
            analysis_id=aid,
            analysis_type=analysis_type,
            project_name=project_name,
            sample_name=sample_names[0] if len(sample_names) == 1 else None,
        )
        # Set up dependencies
        first_j.depends_on(sm_in_progress_j)
        sm_completed_j.depends_on(last_j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        last_j = sm_completed_j
        return last_j
    

def _parse_analysis(data: Dict) -> Optional[Analysis]:
    if not data:
        return None
    if 'id' not in data:
        logger.error(f'Analysis data doesn\'t have id: {data}')
        return None
    if 'type' not in data:
        logger.error(f'Analysis data doesn\'t have type: {data}')
        return None
    if 'status' not in data:
        logger.error(f'Analysis data doesn\'t have status: {data}')
        return None
    a = Analysis(
        id=data['id'],
        type=data['type'],
        status=data['status'],
        sample_ids=set(data.get('sample_ids', [])),
        output=data.get('output', None),
    )
    return a


def parse_reads_from_metadata(  # pylint: disable=too-many-return-statements
    meta: Dict,
    check_existence: bool = True,
) -> Optional[AlignmentInput]:
    """
    Verify the meta.reads object in a sequence db entry
    """
    reads_data = meta.get('reads')
    reads_type = meta.get('reads_type')
    if not reads_data:
        logger.error(f'No "meta/reads" field in meta')
        return None
    if not reads_type:
        logger.error(f'No "meta/reads_type" field in meta')
        return None
    supported_types = ('fastq', 'bam', 'cram')
    if reads_type not in supported_types:
        logger.error(f'ERROR: "reads_type" is expected to be one of {supported_types}')
        return None

    if reads_type in ('bam', 'cram'):
        if len(reads_data) > 1:
            logger.error('Supporting only single bam/cram input')
            return None

        bam_path = reads_data[0]['location']
        if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
            logger.error(
                f'ERROR: expected the file to have an extention .cram or .bam,'
                f'got: {bam_path}'
            )
            return None
        if check_existence and not utils.file_exists(bam_path):
            logger.error(f'ERROR: index file doesn\'t exist: {bam_path}')
            return None

        # Index:
        if not reads_data[0].get('secondaryFiles'):
            logger.error(
                f'ERROR: bam/cram input is expected to have '
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
                f'ERROR: expected the index file to have an extention '
                f'.crai or .bai, got: {index_path}'
            )
        if check_existence and not utils.file_exists(index_path):
            logger.error(f'ERROR: index file doesn\'t exist: {index_path}')
            return None

        return AlignmentInput(bam_or_cram_path=bam_path, index_path=index_path)

    else:
        fqs1 = []
        fqs2 = []
        for lane_data in reads_data:
            assert len(lane_data) == 2, lane_data
            if check_existence and not utils.file_exists(lane_data[0]['location']):
                logger.error(
                    f'ERROR: read 1 file doesn\'t exist: {lane_data[0]["location"]}'
                )
                return None
            if check_existence and not utils.file_exists(lane_data[1]['location']):
                logger.error(
                    f'ERROR: read 2 file doesn\'t exist: {lane_data[1]["location"]}'
                )
                return None

            fqs1.append(lane_data[0]['location'])
            fqs2.append(lane_data[1]['location'])
        return AlignmentInput(fqs1=fqs1, fqs2=fqs2)
