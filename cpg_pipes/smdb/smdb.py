"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import traceback
from textwrap import dedent
from typing import List, Dict, Optional, Collection

from hailtop.batch import Batch
from hailtop.batch.job import Job

from sample_metadata import models
from sample_metadata.apis import (
    SampleApi,
    SequenceApi,
    AnalysisApi,
    ParticipantApi,
)
from sample_metadata.exceptions import ApiException

from cpg_pipes import buckets, images
from cpg_pipes.namespace import Namespace
from cpg_pipes.smdb.types import Analysis, SmSequence

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class SMDB:
    """
    Abstracting the communication with the SampleMetadata database.
    """

    def __init__(
        self, 
        analysis_project: str, 
        do_update_analyses: bool = True,
        do_check_seq_existence: bool = True,
    ):
        """
        :param analysis_project: project where to create the "analysis" entries.
        :param do_update_analyses: if not set, won't update "analysis" entries' 
            statuses.
        :param do_check_seq_existence: when querying "sequence" or "analysis" entries
            with files, check those files existence with gsutil. For "sequence", will
            throw an error. For "analysis", will invalidate by setting status=failure.
        """
        self.sapi = SampleApi()
        self.aapi = AnalysisApi()
        self.seqapi = SequenceApi()
        self.seqapi = SequenceApi()
        self.papi = ParticipantApi()
        self.analysis_project = analysis_project
        self.do_update_analyses = do_update_analyses
        self.do_check_seq_existence = do_check_seq_existence

    def get_samples_by_project(
        self,
        project_names: List[str],
        namespace: Namespace,
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
    ) -> Dict[str, List[Dict]]:
        """
        Returns a dictionary of samples per input projects
        """
        samples_by_project: Dict[str, List[Dict]] = dict()
        for proj_name in project_names:
            input_proj_name = proj_name
            if namespace != Namespace.MAIN:
                input_proj_name += '-test'
            logger.info(f'Finding samples for project {input_proj_name}...')
            samples = self.sapi.get_samples(
                body_get_samples_by_criteria_api_v1_sample_post={
                    'project_ids': [input_proj_name],
                    'active': True,
                }
            )
            logger.info(f'Finding samples for project {input_proj_name}: found {len(samples)}')
            
            participant_id_by_cpgid = self._get_participant_id_by_sid(proj_name)
            
            samples_by_project[proj_name] = []
            for s in samples:
                s['id'] = s['id'].strip()
                s['external_id'] = s['external_id'].strip()
                
                if only_samples:
                    if s['id'] in only_samples or s['external_id'] in only_samples:
                        logger.info(f'Taking sample: {s["id"]}')
                    else:
                        continue
                if skip_samples:
                    if s['id'] in skip_samples or s['external_id'] in skip_samples:
                        logger.info(f'Skiping sample: {s["id"]}')
                        continue
                samples_by_project[proj_name].append(s)
                
                s['participant_id'] = participant_id_by_cpgid.get(s['id'])
                
        return samples_by_project

    def _get_participant_id_by_sid(self, proj_name: str) -> Dict[str, str]:
        try:
            pid_sid = self.papi.get_external_participant_id_to_internal_sample_id(
                proj_name
            )
        except ApiException:
            participant_id_by_cpgid = {}
        else:
            participant_id_by_cpgid = {sid.strip(): pid.strip() for pid, sid in pid_sid}
        return participant_id_by_cpgid

    def find_seq_by_sid(self, sample_ids) -> Dict[str, SmSequence]:
        """
        Return a dict of "Sequence" entries by sample ID
        """
        try:
            seq_infos: List[Dict] = self.seqapi.get_sequences_by_sample_ids(sample_ids)
        except ApiException:
            traceback.print_exc()
            return {}
        else:
            seqs = [SmSequence.parse(d, self) for d in seq_infos]
            seqs_by_sid = {seq.sample_id: seq for seq in seqs}
            return seqs_by_sid

    def update_analysis(self, analysis: 'Analysis', status: str):
        """
        Update "status" of an Analysis entry
        """
        if not self.do_update_analyses:
            return
        try:
            self.aapi.update_analysis_status(
                analysis.id, models.AnalysisUpdateModel(
                    status=models.AnalysisStatus(status)
                )
            )
        except ApiException:
            traceback.print_exc()
        analysis.status = status

    def find_joint_calling_analysis(
        self,
        sample_ids: Collection[str],
    ) -> Optional['Analysis']:
        """
        Query the DB to find the last completed joint-calling analysis for the samples
        """
        try:
            data = self.aapi.get_latest_complete_analysis_for_type(
                project=self.analysis_project,
                analysis_type=models.AnalysisType('joint-calling'),
            )
        except ApiException:
            return None
        a = Analysis.parse(data)
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
        analysis_status: str = 'completed',
        project: Optional[str] = None,
        meta: Optional[Dict] = None
    ) -> Dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type and samples,
        one Analysis object per sample. Assumes the analysis is defined for a single
        sample (e.g. cram, gvcf)
        """
        project = project or self.analysis_project
        analysis_per_sid: Dict[str, Analysis] = dict()

        logger.info(
            f'Querying {analysis_type} analysis entries for project {project}...'
        )
        datas = self.aapi.query_analyses(
            models.AnalysisQueryModel(
                projects=[project],
                sample_ids=sample_ids,
                type=models.AnalysisType(analysis_type),
                status=models.AnalysisStatus(analysis_status),
                meta=meta or {},
            )
        )

        for data in datas:
            a = Analysis.parse(data)
            if not a:
                continue
            assert a.status == 'completed', data
            assert a.type == analysis_type, data
            assert len(a.sample_ids) == 1, data
            analysis_per_sid[list(a.sample_ids)[0]] = a
        logger.info(
            f'Querying {analysis_type} analysis entries for project {project}: found {len(analysis_per_sid)}'
        )
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
        analysis_id: int,
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
        j.image(images.SM_IMAGE)
        cmd = dedent(f"""\
        export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
        gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

        export SM_DEV_DB_PROJECT={self.analysis_project}
        export SM_ENVIRONMENT=PRODUCTION
        
        cat <<EOT >> update.py
        from sample_metadata.apis import AnalysisApi
        from sample_metadata.models import AnalysisUpdateModel, AnalysisStatus
        from sample_metadata import exceptions
        import traceback
        aapi = AnalysisApi()
        try:
            aapi.update_analysis_status(
                analysis_id={analysis_id},
                analysis_update_model=AnalysisUpdateModel(
                    status=AnalysisStatus('{status}')
                ),
            )
        except exceptions.ApiException:
            traceback.print_exc()
        EOT
        python update.py
        """)
        j.command(cmd)
        return j

    def create_analysis(
        self,
        output: str,
        type_: str,
        status: str,
        sample_ids: Collection[str],
        project_name: Optional[str] = None,
    ) -> Optional[int]:
        """
        Tries to create an Analysis entry, returns its id if successfuly
        """
        if not self.do_update_analyses:
            return None

        am = models.AnalysisModel(
            type=models.AnalysisType(type_),
            status=models.AnalysisStatus(status),
            output=output,
            sample_ids=list(sample_ids),
        )
        try:
            aid = self.aapi.create_new_analysis(
                project=project_name or self.analysis_project, analysis_model=am
            )
        except ApiException:
            traceback.print_exc()
            return None
        else:
            logger.info(f'Created analysis of type={type_}, status={status} with ID: {aid}')
            return aid

    def process_existing_analysis(
        self,
        sample_ids: Collection[str],
        completed_analysis: Optional['Analysis'],
        analysis_type: str,
        expected_output_fpath: str,
        project_name: Optional[str] = None,
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
        :param project_name: project name to create new analysis in
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
            elif not buckets.file_exists(found_output_fpath):
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
        if buckets.file_exists(expected_output_fpath):
            logger.info(
                f'Output file {expected_output_fpath} already exists, so creating '
                f'an analysis {label} with status=completed'
            )
            self.create_analysis(
                type_=analysis_type,
                output=expected_output_fpath,
                status='completed',
                sample_ids=sample_ids,
                project_name=project_name,
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
        first_j,
        last_j,
        depends_on,
        project_name: Optional[str] = None,
    ) -> Job:
        if not self.do_update_analyses:
            return last_j
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = self.create_analysis(
            type_=analysis_type,
            output=output_path,
            status='queued',
            sample_ids=sample_names,
            project_name=project_name,
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
        for fj in (first_j if isinstance(first_j, list) else [first_j]):
            fj.depends_on(sm_in_progress_j)
        sm_completed_j.depends_on(last_j)
        if depends_on:
            sm_in_progress_j.depends_on(*depends_on)
        last_j = sm_completed_j
        return last_j
