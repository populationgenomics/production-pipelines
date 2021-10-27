from typing import Collection, List, Dict, Optional, Set
import logging
import traceback

from sample_metadata import (
    AnalysisApi,
    SequenceApi,
    SampleApi,
    AnalysisUpdateModel,
    AnalysisModel,
    exceptions,
)


logger = logging.getLogger(__file__)
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger.setLevel(logging.INFO)


class SMDB:
    """
    Singleton class abstracting the communication with
    the SampleMetadata database
    """

    sapi = SampleApi()
    aapi = AnalysisApi()
    seqapi = SequenceApi()
    do_update_analyses: bool = False

    @classmethod
    def get_samples_by_project(
        cls,
        project: str,
        skip_samples: Optional[Set[str]] = None,
    ) -> List[Dict]:
        """
        Returns a dictionary of samples per input projects
        """
        logger.info(f"Finding samples for project {project}")
        samples = cls.sapi.get_samples(
            body_get_samples_by_criteria_api_v1_sample_post={
                "project_ids": [project],
                "active": True,
            }
        )
        if skip_samples:
            samples = list(filter(lambda s: s["id"] not in skip_samples, samples))

        return samples

    @classmethod
    def get_cram_analyses_from_project_for_sample_ids(
        cls, project: str, sample_ids: List[str]
    ):
        analyses = cls.aapi.get_latest_analysis_for_samples_and_type(
            analysis_type="cram", project=project, request_body=sample_ids
        )
        return analyses

    @classmethod
    def get_variants_for_sample_id_by_sv_type(
        cls, project: str, sample_ids: List[str], sv_type: str
    ):
        # NOTE: this is a proposal, this endpoint does not currently exist like this
        analyses = cls.aapi.get_latest_analysis_for_samples_and_type(
            analysis_type="variants",
            project=project,
            request_body={
                "sample_ids": sample_ids,
                "meta": {
                    # TODO: define some better convention
                    "sv_type": sv_type
                },
            },
        )
        return analyses

    @classmethod
    def create_analysis(
        cls,
        project: str,
        type_: str,
        output: str,
        status: str,
        sample_ids: Collection[str],
            meta: Dict
    ) -> Optional[int]:
        """
        Tries to create an Analysis entry, returns its id if successfuly
        """
        if not cls.do_update_analyses:
            return None

        am = AnalysisModel(
            type=type_,
            output=output,
            status=status,
            sample_ids=sample_ids,
            meta=meta
        )
        aid = cls.aapi.create_new_analysis(project=project, analysis_model=am)
        logger.info(f'Created analysis of type={type_} status={status} with ID: {aid}')
        return aid