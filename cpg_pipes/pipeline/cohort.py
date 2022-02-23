"""
Represents a "cohort" target - all samples from all projects in the pipeline
"""

from os.path import join
from typing import List, Optional
import logging

from cpg_pipes import buckets
from cpg_pipes.hb.inputs import AlignmentInput
from cpg_pipes.pipeline.sample import Sample, PedigreeInfo, Sex
from cpg_pipes.pipeline.project import Project
from cpg_pipes.pipeline.target import Target
from cpg_pipes.smdb.smdb import SMDB
from cpg_pipes.smdb.types import AnalysisType

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class Cohort(Target):
    """
    Represents a "cohort" target - all samples from all projects in the pipeline
    """

    def __init__(self, name: str, pipeline):
        super().__init__()
        self.name = name
        self.pipeline = pipeline
        self._projects: List[Project] = []

    def __repr__(self):
        return self.name

    @property
    def unique_id(self) -> str:
        """
        Unique ID of cohort as a stage "Target" - can be anything, because
        the cohort is expected to be only one
        """
        return self.name

    def get_projects(self, only_active: bool = True) -> List[Project]:
        """
        Gets list of all projects.
        Include only "active" projects (unless only_active is False)
        """
        return [p for p in self._projects if (p.active or not only_active)]

    def get_all_samples(self, only_active: bool = True) -> List[Sample]:
        """
        Gets a flat list of all samples from all projects.
        Include only "active" samples (unless only_active is False)
        """
        all_samples = []
        for proj in self.get_projects(only_active=only_active):
            all_samples.extend(proj.get_samples(only_active))
        return all_samples

    def get_all_sample_ids(self) -> List[str]:
        """
        Gets a flat list of CPG IDs for all samples from all projects.
        """
        return [s.id for s in self.get_all_samples()]

    def add_project(self, name: str) -> Project:
        project_by_name = {p.name: p for p in self._projects}
        if name in project_by_name:
            logger.warning(f'Project {name} already exists')
            return project_by_name[name]
        p = Project(pipeline=self.pipeline, name=name)
        self._projects.append(p)
        return p

    def populate(
        self,
        smdb: SMDB,
        input_projects: List[str],
        local_tmp_dir: str,
        source_tag: Optional[str] = None,
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        ped_files: Optional[List[str]] = None,
        forced_samples: List[str] = None,
    ) -> None:
        """
        Finds input samples, analyses and sequences from the DB,
        populates self.projects, adds pedigree information
        """
        self._populate_projects(
            smdb=smdb,
            input_projects=input_projects,
            skip_samples=skip_samples,
            only_samples=only_samples,
            forced_samples=forced_samples,
            source_tag=source_tag,
        )
        self._populate_seq(smdb)
        self._populate_analysis(smdb, source_tag=source_tag)
        if ped_files:
            self._populate_pedigree(ped_files, local_tmp_dir)

    def _populate_projects(
        self,
        smdb: SMDB,
        input_projects: List[str],
        skip_samples: Optional[List[str]] = None,
        only_samples: Optional[List[str]] = None,
        forced_samples: Optional[List[str]] = None,
        source_tag: Optional[str] = None,
    ):
        samples_by_project = smdb.get_samples_by_project(
            project_names=input_projects,
            namespace=self.pipeline.namespace,
            skip_samples=skip_samples,
            only_samples=only_samples,
        )
        for proj_name, sample_datas in samples_by_project.items():
            project = self.add_project(name=proj_name)
            for s_data in sample_datas:
                meta = s_data.get('meta', {})
                if source_tag:
                    meta['source_tag'] = source_tag
                external_id = s_data['external_id']
                participant_id = s_data.get('participant_id')
                s = project.add_sample(
                    id=s_data['id'],
                    external_id=external_id.strip(),
                    participant_id=participant_id.strip() if participant_id else None,
                    **s_data.get('meta', dict()),
                )
                if forced_samples and s.id in forced_samples:
                    logger.info(f'Force rerunning sample {s.id} even if outputs exist')
                    s.forced = True

    def _populate_pedigree(self, ped_files: List[str], local_tmp_dir: str):
        sample_by_participant_id = dict()
        for s in self.get_all_samples():
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
        for project in self.get_projects():
            samples_with_ped = [s for s in project.get_samples() if s.pedigree]
            logger.info(
                f'{project.name}: found pedigree info for {len(samples_with_ped)} '
                f'samples out of {len(project.get_samples())}'
            )

    def _populate_seq(self, smdb: SMDB):
        """
        Queries Sequence entries for each sample
        """
        all_sample_ids = self.get_all_sample_ids()
        seqs_by_sid = smdb.find_seq_by_sid(all_sample_ids)
        for s in self.get_all_samples():
            if s.id in seqs_by_sid:
                s.seq = seqs_by_sid[s.id]
                assert s.seq is not None
                s.alignment_input = s.seq.parse_reads()

    def _populate_analysis(self, smdb: SMDB, source_tag: Optional[str] = None):
        all_sample_ids = self.get_all_sample_ids()

        jc_analysis = smdb.find_joint_calling_analysis(
            sample_ids=all_sample_ids,
        )

        for project in self.get_projects():
            sample_ids = [s.id for s in project.get_samples()]

            cram_per_sid = smdb.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type=AnalysisType.CRAM.value,
                meta={'source': source_tag} if source_tag else None,
                project=project.name,
            )
            gvcf_per_sid = smdb.find_analyses_by_sid(
                sample_ids=sample_ids,
                analysis_type=AnalysisType.GVCF.value,
                meta={'source': source_tag} if source_tag else None,
                project=project.name,
            )

            for s in project.get_samples():
                if s.id in cram_per_sid:
                    s.analysis_by_type[AnalysisType.CRAM] = cram_per_sid[s.id]
                    cram_path = cram_per_sid[s.id].output
                    if cram_path:
                        index_path = cram_path + '.crai'
                        s.alignment_input = AlignmentInput(
                            bam_or_cram_path=cram_path, index_path=index_path
                        )
                if s.id in gvcf_per_sid:
                    s.analysis_by_type[AnalysisType.GVCF] = gvcf_per_sid[s.id]

        if jc_analysis:
            self.analysis_by_type[AnalysisType.JOINT_CALLING] = jc_analysis
