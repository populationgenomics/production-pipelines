"""
Helpers to pull intermediate results from previous Batch runs.
"""

import logging
from dataclasses import dataclass
from typing import Optional, Tuple, Dict

from cloudpathlib import CloudPath

from cpg_pipes import buckets
from cpg_pipes.hb.batch import get_hail_bucket

logger = logging.getLogger(__file__)


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
        tmp_bucket: CloudPath,
        keep_scratch: bool,
    ) -> Dict[Tuple[Optional[str], str], 'PrevJob']:
        """
        `fpath` is the path to result.txt, generated from mysql (see above), 
        and `batchid` is a 6-letter batch ID (e.g. feb0e9 
        in gs://cpg-seqr-main-tmp/seqr_loader/v0/hail/batch/feb0e9)
        which is not to be confused with a batch number from the URL.

        Returns a dictionary of PrevJob objects indexed by a pair of optional 
        CPG ID string (for sample-level jobs) and a job name string.
        """
        hail_bucket = get_hail_bucket(tmp_bucket, keep_scratch)
        
        if not buckets.exists(fpath):
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
