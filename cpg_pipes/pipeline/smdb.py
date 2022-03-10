"""
Functions to find the pipeline inputs and communicate with the SM server
"""

import logging
import traceback
from os.path import join
from textwrap import dedent

from cloudpathlib import CloudPath
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

from .. import buckets, images
from .cohort import Cohort
from .sample import PedigreeInfo, Sex, Sample
from .dataset import Dataset
from .analysis import Analysis, AnalysisType, AnalysisStatus
from .sequence import SmSequence

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class SMDB:
    """
    Communication with the SampleMetadata database.
    """

    def __init__(
        self, 
        analysis_dataset: str, 
        do_update_analyses: bool = True,
        do_check_seq_existence: bool = True,
    ):
        """
        :param analysis_dataset: dataset where to create the "analysis" entries.
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
        self.analysis_dataset = analysis_dataset
        self.do_update_analyses = do_update_analyses
        self.do_check_seq_existence = do_check_seq_existence

    def populate_cohort(
        self,
        cohort: Cohort,
        local_tmp_dir: str,
        skip_samples: list[str]|None = None,
        only_samples: list[str]|None = None,
        ped_files: list[str]|None = None,
    ) -> Cohort:
        """
        Populate database entries for all datasets.
        """
        for dataset in cohort.get_datasets():
            self.populate_dataset(
                dataset,
                only_samples=only_samples,
                skip_samples=skip_samples,
            )
        self._populate_analysis(cohort)
        if ped_files:
            self._populate_pedigree(cohort, ped_files, local_tmp_dir)
        return cohort
    
    def _populate_pedigree(
        self, 
        cohort: Cohort, 
        ped_files: list[str], 
        local_tmp_dir: str
    ):
        sample_by_participant_id = dict()
        for s in cohort.get_all_samples():
            sample_by_participant_id[s.participant_id] = s

        for i, ped_file in enumerate(ped_files):
            local_ped_file = join(local_tmp_dir, f'ped_file_{i}.ped')
            buckets.gsutil_cp(ped_file, local_ped_file)
            with open(local_ped_file) as f:
                for line in f:
                    fields = line.strip().split('\t')[:6]
                    fam_id, sam_id, pat_id, mat_id, sex, phenotype = fields
                    if sam_id in sample_by_participant_id:
                        s = sample_by_participant_id[sam_id]
                        s.pedigree = PedigreeInfo(
                            sample=s,
                            fam_id=fam_id,
                            dad=sample_by_participant_id.get(pat_id),
                            mom=sample_by_participant_id.get(mat_id),
                            sex={
                                '1': Sex.MALE, 
                                '2': Sex.FEMALE,
                                'M': Sex.MALE,
                                'F': Sex.FEMALE,
                            }.get(sex, Sex.UNKNOWN),
                            phenotype=phenotype or '0',
                        )
        for dataset in cohort.get_datasets():
            samples_with_ped = [s for s in dataset.get_samples() if s.pedigree]
            logger.info(
                f'{dataset.name}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(dataset.get_samples())}'
            )

    def _populate_analysis(self, cohort: Cohort, source_tag: str|None = None):
        all_sample_ids = cohort.get_all_sample_ids()

        jc_analysis = self.find_joint_calling_analysis(
            sample_ids=all_sample_ids,
        )
        if jc_analysis:
            cohort.analysis_by_type[AnalysisType.JOINT_CALLING] = jc_analysis

        for dataset in cohort.get_datasets():
            sample_ids = [s.id for s in dataset.get_samples()]

            for atype in AnalysisType:
                analysis_per_sid = self.find_analyses_by_sid(
                    sample_ids=sample_ids,
                    analysis_type=atype,
                    meta={'source': source_tag} if source_tag else None,
                    dataset=dataset.name,
                )
                for s in dataset.get_samples():
                    if s.id in analysis_per_sid:
                        s.analysis_by_type[atype] = analysis_per_sid[s.id]

    def populate_dataset(
        self,
        dataset: Dataset,
        skip_samples: list[str]|None = None,
        only_samples: list[str]|None = None,
    ) -> Dataset:
        """
        Populate database entries for one dataset.
        """
        sample_datas = self.__get_sample_entries(
            dataset_name=dataset.name,
            skip_samples=skip_samples,
            only_samples=only_samples,
        )

        participant_id_by_cpgid = self.__get_participant_id_by_sid(dataset.name)

        samples = [
            dataset.add_sample(
                id=sample_d['id'], 
                external_id=sample_d['external_id'],
                participant_id=participant_id_by_cpgid.get(sample_d['id']),
                **sample_d.get('meta', {}),
            )
            for sample_d in sample_datas
        ]

        self.populate_sequence(samples)

        return dataset

    def populate_sequence(self, samples: list[Sample]) -> list[Sample]:
        """
        Get and parse "Sequence" entries for each sample.
        """
        try:
            seq_infos: list[dict] = self.seqapi.get_sequences_by_sample_ids(
                [s.id for s in samples]
            )
        except ApiException:
            traceback.print_exc()
            return []
        
        seq_info_by_sid = {seq['sample_id']: seq for seq in seq_infos}
        for sample in samples:
            seq_info = seq_info_by_sid[sample.id]
            sample.seq = SmSequence.parse(seq_info, self.do_check_seq_existence) 
        return samples
    
    def __get_sample_entries(
        self,
        dataset_name: str,
        skip_samples: list[str]|None = None,
        only_samples: list[str]|None = None,
    ) -> list[dict]:
        """
        Get Sample entries and apply `skip_samples` and `only_samples` filters.
        This is a helper method; use public `populate_dataset` to parse samples. 
        """
        logger.info(f'Finding samples for dataset {dataset_name}...')
        sample_entries = self.sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                'project_ids': [dataset_name],
                'active': True,
            }
        )
        logger.info(
            f'Finding samples for dataset {dataset_name}:'
            f'found {len(sample_entries)}'
        )
        
        filtered_entries = []
        for sample_dict in sample_entries:
            cpgid = sample_dict['id'].strip()
            extid = sample_dict['external_id'].strip()

            if only_samples:
                if cpgid in only_samples or extid in only_samples:
                    logger.info(f'Picking sample: {cpgid}|{extid}')
                else:
                    continue
            if skip_samples:
                if cpgid in skip_samples or extid in skip_samples:
                    logger.info(f'Skiping sample: {cpgid}|{extid}')
                    continue 
            filtered_entries.append(sample_dict)
        return filtered_entries

    def __get_participant_id_by_sid(self, proj_name: str) -> dict[str, str]:
        try:
            pid_sid = self.papi.get_external_participant_id_to_internal_sample_id(
                proj_name
            )
        except ApiException:
            participant_id_by_cpgid = {}
        else:
            participant_id_by_cpgid = {sid.strip(): pid.strip() for pid, sid in pid_sid}
        return participant_id_by_cpgid

    def update_analysis(self, analysis: 'Analysis', status: AnalysisStatus):
        """
        Update "status" of an Analysis entry
        """
        if not self.do_update_analyses:
            return
        try:
            self.aapi.update_analysis_status(
                analysis.id, models.AnalysisUpdateModel(
                    status=models.AnalysisStatus(status.value)
                )
            )
        except ApiException:
            traceback.print_exc()
        analysis.status = status

    def find_joint_calling_analysis(
        self,
        sample_ids: list[str],
    ) -> Analysis|None:
        """
        Query the DB to find the last completed joint-calling analysis for the samples
        """
        try:
            data = self.aapi.get_latest_complete_analysis_for_type(
                project=self.analysis_dataset,
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
        dataset: str|None = None,
        meta: dict|None = None
    ) -> dict[str, Analysis]:
        """
        Query the DB to find the last completed analysis for the type and samples,
        one Analysis object per sample. Assumes the analysis is defined for a single
        sample (e.g. cram, gvcf)
        """
        dataset = dataset or self.analysis_dataset
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
            f'Querying {analysis_type} analysis entries for dataset {dataset}: found {len(analysis_per_sid)}'
        )
        return analysis_per_sid

    def make_sm_in_progress_job(self, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "in-progress"
        """
        kwargs['status'] = AnalysisStatus.IN_PROGRESS.value
        return self.make_sm_update_status_job(*args, **kwargs)

    def make_sm_completed_job(self, *args, **kwargs) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status
        to "completed"
        """
        kwargs['status'] = AnalysisStatus.COMPLETED.value
        return self.make_sm_update_status_job(*args, **kwargs)

    def make_sm_update_status_job(
        self,
        b: Batch,
        analysis_id: int,
        analysis_type: str,
        status: str,
        sample_name: str|None = None,
        dataset_name: str|None = None,
    ) -> Job:
        """
        Creates a job that updates the sample metadata server entry analysis status.
        """
        assert status in ['in-progress', 'failed', 'completed', 'queued']
        job_name = ''
        if dataset_name and sample_name:
            job_name += f'{dataset_name}/{sample_name}: '
        job_name += f'Update SM: {analysis_type} to {status}'

        if not self.do_update_analyses:
            return b.new_job(f'{job_name} [skip]')

        j = b.new_job(job_name)
        j.image(images.SM_IMAGE)
        cmd = dedent(f"""\
        export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
        gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

        export SM_DEV_DB_PROJECT={self.analysis_dataset}
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
        output: CloudPath | str,
        type_: str | AnalysisType,
        status: str | AnalysisStatus,
        sample_ids: list[str],
        dataset_name: str|None = None,
    ) -> int|None:
        """
        Tries to create an Analysis entry, returns its id if successfuly
        """
        if not self.do_update_analyses:
            return None

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
            aid = self.aapi.create_new_analysis(
                project=dataset_name or self.analysis_dataset, analysis_model=am
            )
        except ApiException:
            traceback.print_exc()
            return None
        else:
            logger.info(f'Created analysis of type={type_}, status={status} with ID: {aid}')
            return aid

    def process_existing_analysis(
        self,
        sample_ids: list[str],
        completed_analysis: Analysis | None,
        analysis_type: str,
        expected_output_fpath: CloudPath,
        dataset_name: str | None = None,
    ) -> CloudPath | None:
        """
        Checks whether existing analysis exists, and output matches the expected output
        file. Invalidates bad analysis by setting status=failure, and submits a
        status=completed analysis if the expected output already exists.

        Returns the path to the output if it can be reused, otherwise None.

        :param sample_ids: sample IDs to pull the analysis for
        :param completed_analysis: existing completed analysis of this type for these 
        samples
        :param analysis_type: cram, gvcf, joint_calling
        :param expected_output_fpath: where the pipeline expects the analysis output 
        file to sit on the bucket (will invalidate the analysis if it doesn't match)
        :param dataset_name: the name of the dataset where to create a new analysis
        :return: path to the output if it can be reused, otherwise None
        """
        label = f'type={analysis_type}'
        if len(sample_ids) > 1:
            label += f' for {", ".join(sample_ids)}'

        found_output_fpath: CloudPath | None = None
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
            self.update_analysis(completed_analysis, status=AnalysisStatus.FAILED)

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
                dataset_name=dataset_name,
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
        b: Batch,
        analysis_type: AnalysisType,
        output_path: CloudPath | str,
        sample_names: list[str],
        first_j: Job,
        last_j: Job,
        depends_on: list[Job] | None,
        dataset_name: str | None = None,
    ) -> Job:
        """
        Given a first job and a last job object, insert corresponding "queued", 
        "in_progress", and "completed" jobs. Useful to call in the end of a stage
        definition.
        """
        if not self.do_update_analyses:
            return last_j
        # Interacting with the sample metadata server:
        # 1. Create a "queued" analysis
        aid = self.create_analysis(
            type_=analysis_type,
            output=output_path,
            status='queued',
            sample_ids=sample_names,
            dataset_name=dataset_name,
        )
        # 2. Queue a job that updates the status to "in-progress"
        sm_in_progress_j = self.make_sm_in_progress_job(
            b,
            analysis_id=aid,
            analysis_type=analysis_type,
            dataset_name=dataset_name,
            sample_name=sample_names[0] if len(sample_names) == 1 else None,
        )
        # 2. Queue a job that updates the status to "completed"
        sm_completed_j = self.make_sm_completed_job(
            b,
            analysis_id=aid,
            analysis_type=analysis_type,
            dataset_name=dataset_name,
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
