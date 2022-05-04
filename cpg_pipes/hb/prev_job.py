"""
Helpers to pull intermediate results from previous Batch runs.
"""

import logging
from dataclasses import dataclass
from typing import Optional, Tuple, Dict

from .. import Path
from .. import utils

logger = logging.getLogger(__file__)


@dataclass
class PrevJob:
    """
    Represents a Job from a previous Batch, parsed from the Batch database with:
    BATCH=6553
    mysql --ssl-ca=/sql-config/server-ca.pem --ssl-cert=/sql-config/client-cert.pem --ssl-key=/sql-config/client-key.pem --host=10.125.0.3 --user=root --password batch -e "SELECT j.batch_id, j.job_id, ja.key, ja.value FROM jobs as j INNER JOIN job_attributes as ja ON j.job_id = ja.job_id WHERE j.batch_id=$BATCH AND ja.batch_id=j.batch_id AND ja.key='name' AND j.state='Success';" > result.txt
    """

    cpgid: str | None
    projname: str | None
    batch_number: int
    job_number: int
    jtype: str
    jname: str
    batchid: str
    hail_bucket: Path

    @staticmethod
    def parse(
        fpath: Path,
        batchid: str,
        hail_bucket: Path,
    ) -> Dict[Tuple[Optional[str], str], 'PrevJob']:
        """
        `fpath` is the path to result.txt, generated from mysql (see above),
        and `batchid` is a 6-letter batch ID (e.g. feb0e9
        in gs://cpg-seqr-main-tmp/seqr_loader/v0/hail/batch/feb0e9)
        which is not to be confused with a batch number from the URL.

        Returns a dictionary of PrevJob objects indexed by a pair of optional
        CPG ID string (for sample-level jobs) and a job name string.
        """
        if not utils.exists(fpath):
            return dict()
        prev_batch_jobs: Dict[Tuple[Optional[str], str], PrevJob] = dict()
        with fpath.open() as f:
            for line in f:
                if 'batch_id' in line.split():
                    continue  # skip header
                batch_number, job_number, _, jname = line.strip().split('\t')

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
