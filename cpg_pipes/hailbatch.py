from dataclasses import dataclass
from typing import Optional, Union, List, Tuple, Dict
import hailtop.batch as hb

from cpg_pipes import utils


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
    cpgid: str
    projname: str
    batch_number: int
    job_number: int
    jtype: str
    jname: str
    batchid: int
    hail_bucket: str

    @staticmethod
    def parse(
        fpath: str, 
        batchid: str, 
        hail_bucket: str,
    ) -> Dict:
        """
        fpath is the path to result.txt (see above), and batchid is a 6-letter
        batch ID, e.g. feb0e9 in gs://cpg-seqr-main-tmp/seqr_loader/v0/hail/batch/feb0e9
        which is not to be confused with a batch number from the URL
        """
        if not utils.file_exists(fpath):
            return dict()
        prev_batch_jobs = dict()
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
                    batch_number=batch_number,
                    job_number=job_number,
                    jtype=jtype,
                    jname=jname,
                    batchid=batchid,
                    hail_bucket=hail_bucket,
                )
        return prev_batch_jobs
