import logging
import os
from dataclasses import dataclass
from typing import Optional, Union, List, Tuple, Dict
import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_pipes import utils

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)



@dataclass
class BamOrCramAlignmentInput:
    """
    Represents inputs for an alignment job of CRAM or BAM files
    """
    bam_or_cram_path: Optional[Union[str, hb.ResourceGroup]] = None
    index_path: Optional[str] = None

    def get_cram_input_group(self, b) -> hb.ResourceGroup:
        """
        Makes a ResourceGroup of bam/cram with accompanying index
        """
        assert self.bam_or_cram_path
        if isinstance(self.bam_or_cram_path, str):
            return b.read_input_group(
                base=self.bam_or_cram_path,
                index=self.index_path or self.bam_or_cram_path + '.crai'
            )
        else:
            return self.bam_or_cram_path


@dataclass
class FqAlignmentInput:
    """
    Represents inputs for an alignment job of fastq files
    """
    fqs1: Optional[List[Union[str, hb.ResourceFile]]] = None
    fqs2: Optional[List[Union[str, hb.ResourceFile]]] = None

    def get_fq_inputs(self, b) -> Tuple[List[hb.Resource], List[hb.Resource]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        assert self.fqs1 is not None and self.fqs2 is not None 
        if isinstance(self.fqs1[0], hb.Resource):
            files1 = self.fqs1
            files2 = self.fqs2
        else:
            files1 = [b.read_input(f1) for f1 in self.fqs1]
            files2 = [b.read_input(f1) for f1 in self.fqs2]
        return files1, files2


@dataclass
class AlignmentInput(BamOrCramAlignmentInput, FqAlignmentInput):
    pass


@dataclass
class PrevJob:
    """
    Represents a Job from a previous Batch, parsed from the Batch database with:
    BATCH=6553
    mysql --ssl-ca=/sql-config/server-ca.pem --ssl-cert=/sql-config/client-cert.pem --ssl-key=/sql-config/client-key.pem --host=10.125.0.3 --user=root --password batch -e "SELECT j.batch_id, j.job_id, ja.key, ja.value FROM jobs as j INNER JOIN job_attributes as ja ON j.job_id = ja.job_id WHERE j.batch_id=$BATCH AND ja.batch_id=j.batch_id AND ja.key='name' AND j.state='Success';" > result.txt    
    """
    cpgid: Optional[str]
    projname: Optional[str]
    batch_number: int
    job_number: int
    jtype: str
    jname: str
    batchid: str
    hail_bucket: str

    @staticmethod
    def parse(
        fpath: str, 
        batchid: str, 
        hail_bucket: str,
    ) -> Dict[Tuple[Optional[str], str], 'PrevJob']:
        """
        `fpath` is the path to result.txt, generated from mysql (see above), 
        and `batchid` is a 6-letter batch ID (e.g. feb0e9 
        in gs://cpg-seqr-main-tmp/seqr_loader/v0/hail/batch/feb0e9)
        which is not to be confused with a batch number from the URL.

        Returns a dictionary of PrevJob objects indexed by a pair of optional 
        CPG ID string (for sample-level jobs) and a job name string.
        """
        if not utils.file_exists(fpath):
            return dict()
        prev_batch_jobs: Dict[Tuple[Optional[str], str], PrevJob] = dict()
        with open(fpath) as f:
            for line in f:
                if 'batch_id' in line.split():
                    continue  # skip header
                batch_number, job_number, _, jname = \
                    line.strip().split('\t')

                proj_cpgid = ''
                if ': ' in jname:
                    proj_cpgid, jname = jname.split(': ')
                jtype = jname.split()[0]
                projname, cpgid = None, None
                if '/' in proj_cpgid:
                    projname, cpgid = proj_cpgid.split('/')

                key = (cpgid, jname)
                assert key not in prev_batch_jobs, (key, prev_batch_jobs)
                prev_batch_jobs[key] = PrevJob(
                    cpgid=cpgid,
                    projname=projname,
                    batch_number=int(batch_number),
                    job_number=int(job_number),
                    jtype=jtype,
                    jname=jname,
                    batchid=batchid,
                    hail_bucket=hail_bucket,
                )
        return prev_batch_jobs


class Batch(hb.Batch):
    """
    Inheriting from Hail `Batch` class just so we can register added jobs to
    print statistics before submitting.
    """
    def __init__(self, name, backend, *args, **kwargs):
        super().__init__(name, backend, *args, **kwargs)
        # Job stats registry
        self.labelled_jobs = dict()
        self.other_job_num = 0
        self.total_job_num = 0

    def new_job(
        self,
        name: Optional[str] = None,
        attributes: Optional[Dict[str, str]] = None,
        **kwargs,
    ) -> Job:
        """
        Adds job to the Batch, and also registers it in `self.job_stats` for
        statistics.
        """
        if not name:
            logger.critical('Error: job name must be defined')
        
        attributes = attributes or dict()
        project = attributes.get('project')
        sample = attributes.get('sample')
        samples = attributes.get('samples')
        label = attributes.get('label', name)

        if sample and project:
            name = f'{project}/{sample}: {name}'
        elif project:
            name = f'{project}: {name}'

        if label and (sample or samples):
            if label not in self.labelled_jobs:
                self.labelled_jobs[label] = {'job_n': 0, 'samples': set()}
            self.labelled_jobs[label]['job_n'] += 1
            self.labelled_jobs[label]['samples'] |= (samples or {sample})
        else:
            self.other_job_num += 1
        self.total_job_num += 1
        j = super().new_job(name, attributes=attributes)
        return j
    

def get_hail_bucket(tmp_bucket, keep_scratch):
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    return hail_bucket


def setup_batch(
    title: str, 
    keep_scratch: bool,
    tmp_bucket: str,
    analysis_project_name: str,
    billing_project: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> Batch:
    if not hail_bucket:
        hail_bucket = get_hail_bucket(tmp_bucket, keep_scratch)
    billing_project = (
        billing_project or
        os.getenv('HAIL_BILLING_PROJECT') or
        analysis_project_name
    )
    logger.info(
        f'Starting Hail Batch with the project {billing_project}, '
        f'bucket {hail_bucket}'
    )
    backend = hb.ServiceBackend(
        billing_project=billing_project,
        bucket=hail_bucket.replace('gs://', ''),
        token=os.environ['HAIL_TOKEN'],
    )
    b = Batch(name=title, backend=backend)
    return b
