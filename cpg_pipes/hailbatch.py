"""
Utility functions to interact with Hail Batch Jobs and Resources
"""

import logging
import math
import os
from dataclasses import dataclass
from typing import Optional, Union, List, Tuple, Dict, cast
import hailtop.batch as hb
from hailtop.batch import ResourceGroup
from hailtop.batch.job import Job

from cpg_pipes import utils, resources

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


@dataclass
class AlignmentInput:
    """
    Represents inputs for an alignment job, which can be a set of fastq files,
    or a CRAM or a BAM file with an index.
    """
    fqs1: Optional[List[Union[str, hb.ResourceFile]]] = None
    fqs2: Optional[List[Union[str, hb.ResourceFile]]] = None
    bam_or_cram_path: Optional[Union[str, hb.ResourceGroup]] = None
    index_path: Optional[str] = None
    
    def is_fastq(self) -> bool:
        """
        Checks that it's a fastq pair, and both in pair are of the same type and length
        """
        if self.fqs1 or self.fqs2:
            assert self.fqs1 and self.fqs2, self
            if any(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, str) for fq in [self.fqs1, self.fqs2]), self
            elif any(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]):
                assert all(isinstance(fq, hb.ResourceFile) for fq in [self.fqs1, self.fqs2]), self
            else:
                assert len(self.fqs1) == len(self.fqs2), self
            return True
        assert self.bam_or_cram_path, self
        return False

    def is_bam_or_cram(self) -> bool:
        """
        Checks that it's a BAM or a CRAM file
        """
        if self.bam_or_cram_path:
            return True
        assert self.fqs1 and self.fqs2, self
        return False

    def get_fqs1(self) -> List[Union[str, hb.ResourceFile]]:
        assert self.is_fastq()
        return cast(List, self.fqs1)

    def get_fqs2(self) -> List[Union[str, hb.ResourceFile]]:
        assert self.is_fastq()
        return cast(List, self.fqs2)

    def as_fq_inputs(self, b) -> Tuple[List[hb.Resource], List[hb.Resource]]:
        """
        Makes a pair of lists of ResourceFile objects for fqs1 and fqs2
        """
        assert self.is_fastq()
        self.fqs1 = cast(List, self.fqs1)
        self.fqs2 = cast(List, self.fqs2)
        if isinstance(self.fqs1[0], hb.Resource):
            files1 = self.fqs1
            files2 = self.fqs2
        else:
            files1 = [b.read_input(f1) for f1 in self.fqs1]
            files2 = [b.read_input(f1) for f1 in self.fqs2]
        return files1, files2

    def as_cram_input_group(self, b) -> hb.ResourceGroup:
        """
        Makes a ResourceGroup of bam/cram with accompanying index
        """
        assert self.is_bam_or_cram()
        self.bam_or_cram_path = cast(str, self.bam_or_cram_path)
        index_path = self.index_path
        if not index_path:
            if self.bam_or_cram_path.endswith('.bam'):
                index_path = self.bam_or_cram_path + '.bai'
            else:
                assert self.bam_or_cram_path.endswith('.cram'), self.bam_or_cram_path
                index_path = self.bam_or_cram_path + '.crai'
                
        if isinstance(self.bam_or_cram_path, str):
            return b.read_input_group(
                base=self.bam_or_cram_path,
                index=index_path
            )
        else:
            return self.bam_or_cram_path


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
    
    
def job_name(name, sample: str = None, project: str = None):
    if sample and project:
        name = f'{project}/{sample}: {name}'
    elif project:
        name = f'{project}: {name}'
    return name


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

        name = job_name(name, sample, project)

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


def get_hail_bucket(tmp_bucket, keep_scratch) -> str:
    hail_bucket = os.environ.get('HAIL_BUCKET')
    if not hail_bucket or keep_scratch:
        # Scratch files are large, so we want to use the temporary bucket to put them in
        hail_bucket = f'{tmp_bucket}/hail'
    return hail_bucket


def setup_batch(
    title: str, 
    analysis_project_name: str = 'fewgenomes',
    tmp_bucket: str = 'gs://cpg-fewgenomes-test-tmp/hail',
    keep_scratch: bool = True,
    billing_project: Optional[str] = None,
    hail_bucket: Optional[str] = None,
) -> Batch:
    """
    Wrapper around Backend and Batch initialization. 
    Handles setting the tmp bucket and billing project.
    """
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
        token=os.environ.get('HAIL_TOKEN'),
    )
    b = Batch(name=title, backend=backend)
    return b


def fasta_ref_resource(b) -> ResourceGroup:
    """
    Returns fasta reference resource group
    """
    return b.read_input_group(**resources.REF_D)


GCLOUD_CMD = """\
export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account \
--key-file=$GOOGLE_APPLICATION_CREDENTIALS
"""

RETRY_CMD = """\
function fail {
  echo $1 >&2
  exit 1
}

function retry {
  local n_attempts=10
  local delay=30
  local n=1
  while true; do
    "$@" && break || {
      if [[ $n -lt $n_attempts ]]; then
        ((n++))
        echo "Command failed. Attempt $n/$n_attempts:"
        sleep $delay;
      else
        fail "The command has failed after $n attempts."
      fi
    }
  done
}

function retry_gs_cp {
  src=$1

  if [ -n "$2" ]; then
    dst=$2
  else
    dst=/io/batch/${basename $src}
  fi
  
  retry gsutil -o GSUtil:check_hashes=never cp $src $dst
}
"""


def wrap_command(
    command: str,
    monitor_space: bool = False,
    setup_gcp: bool = False,
    define_retry_function: bool = False,
    dedent: bool = True,
) -> str:
    """
    Wraps a command for submission
    If job_resource is defined, monitors output space.
    If output_bucket_path_to_check is defined, checks if this file(s) exists,
    and if it does, skips running the rest of the job.
    
    :param command: command to wrap
    :param monitor_space: add a background process that checks the instance disk 
        space every 5 minutes and prints it to the screen
    :param setup_gcp: login to GCP
    :param define_retry_function: when set, adds bash functions `retry` that attempts 
        to redo a command with a pause of default 30 seconds (useful to pull inputs 
        and get around GoogleEgressBandwidth Quota or other google quotas)
    :param dedent: remove all common leading intendation from the command
    """
    cmd = f"""\
    set -o pipefail
    set -ex
    {GCLOUD_CMD if setup_gcp else ''}
    {RETRY_CMD if define_retry_function else ''}
    
    {f'(while true; do {monitor_space_command()}; sleep 600; done) &'
    if monitor_space else ''}
    
    {command}
    
    {monitor_space_command() if monitor_space else ''}
    """

    if dedent:
        # remove any leading spaces and tabs
        cmd = '\n'.join(line.strip() for line in cmd.split('\n'))
        # remove sretches of spaces
        cmd = '\n'.join(' '.join(line.split()) for line in cmd.split('\n'))
    return cmd


def monitor_space_command():
    """
    Make command that monitors the instance storage space and memory
    """
    return f'df -h; du -sh /io; du -sh /io/batch'


@dataclass(init=False)
class MachineType:
    """
    Hail Batch machine type on GCP
    """
    MIN_NCPU: int = 2
    THREADS_ON_CPU = 2  # hyperthreading

    def __init__(
        self,
        name: str,
        ncpu: int,
        mem_gb_per_core: float,
        price_per_hour: float,
    ):
        self.name = name
        self.max_ncpu = ncpu
        self.mem_gb_per_core = mem_gb_per_core
        self.price_per_hour = price_per_hour
    
    def calc_instance_disk_gb(self) -> int:
        """
        The maximum available storage on an instance is calculated 
        in `batch/batch/utils.py/unreserved_worker_data_disk_size_gib()`
        as the disk size (375G) minus reserved image size (30G) minus
        reserved storage per core (5G*ncpu = 120G for a 32-core instance),
        """
        disk_size_gb = 375
        reserved_gb = 30
        reserved_gb_per_core = 5
        return disk_size_gb - reserved_gb - reserved_gb_per_core * self.max_ncpu

    def set_resources(
        self, 
        j: Job,
        ncpu: Optional[int] = None, 
        nthreads: Optional[int] = None, 
        mem_gb: Optional[float] = None,
        storage_gb: Optional[float] = None, 
        attach_disk_storage_gb: Optional[float] = None,
    ) -> 'JobResource':
        """
        Set resources to a Job object. If any of other parameters are set,
        they will be used as a bound to request a fraction of an instance.
        """
        return self.request_resources(
            ncpu=ncpu,
            nthreads=nthreads,
            mem_gb=mem_gb,
            storage_gb=storage_gb,
            attach_disk_storage_gb=attach_disk_storage_gb,
        ).set_to_job(j)

    def request_resources(
        self,
        ncpu: Optional[int] = None, 
        nthreads: Optional[int] = None, 
        mem_gb: Optional[float] = None,
        storage_gb: Optional[float] = None, 
        attach_disk_storage_gb: Optional[float] = None,
    ) -> 'JobResource':
        """
        Request resources from the machine, satisfying all provided requirements.
        If not requirements are provided, the minimal amount of cores 
        (self.MIN_NCPU) will be used.
        """
        # determining the biggest limit to satisfy, measured in the number of CPUs:
        min_ncpu = max(filter(None, [
            self.adjust_ncpu(ncpu or self.MIN_NCPU),
            self.nthreads_to_ncpu(nthreads) if nthreads else None,
            self.mem_gb_to_ncpu(mem_gb) if mem_gb else None,
            self.storage_gb_to_ncpu(storage_gb) if storage_gb else None,
        ]))
        return JobResource(
            machine_type=self, 
            ncpu=min_ncpu, 
            attach_disk_storage_gb=attach_disk_storage_gb
        )

    def mem_gb_to_ncpu(self, mem_gb: float) -> int:
        """
        Converts memory requirement to the number of CPU requirement.
        """
        ncpu = int(math.ceil(mem_gb / self.mem_gb_per_core))
        return self.adjust_ncpu(ncpu)

    def storage_gb_to_ncpu(self, storage_gb: float) -> int:
        """
        Converts storage requirement to the number of CPU requirement.

        We want to avoid attaching disks (attaching a disk to an existing instance
        might fail with `mkfs.ext4 ...` error, see:
        https://batch.hail.populationgenomics.org.au/batches/7488/jobs/12
        So this function will calculate the number of CPU to request so your jobs 
        can be packed to fit the default instance's available storage 
        (caluclated with self.calc_instance_disk_gb()).
        """
        njobs_fit_on_machine = self.calc_instance_disk_gb() / storage_gb
        return self.adjust_ncpu(int(math.ceil(self.max_ncpu / njobs_fit_on_machine)))

    def nthreads_to_ncpu(self, nthreads: int) -> int:
        return self.adjust_ncpu(math.ceil(nthreads / 2))

    def adjust_ncpu(self, ncpu: int) -> int:
        """
        Adjust request number of CPU to a correct requestable number:
        the nearest power of 2, not less than the minimal number of cores allowed.
        """
        if ncpu > self.max_ncpu:
            ValueError(
                f'Requesting more cores than available on {self.name} machine: '
                f'{ncpu}>{self.max_ncpu}'
            )

        if ncpu < MachineType.MIN_NCPU:
            logger.warning(
                f'align: ncpu is adjusted: {ncpu} -> {MachineType.MIN_NCPU}, '
                f'to the minimal amount to request from an instance.'
            )
            ncpu = MachineType.MIN_NCPU

        # round to the nearest power of 2 (15 -> 16, 16 -> 16, 17 -> 32)
        return int(pow(2, math.ceil(math.log2(ncpu))))
    

# Default Hail Batch machine. 
# 
# Bigger default number of cores (32 vs 16) would allow for more threads for alignment,
# however would mean a smaller available storage (5G per core is reserved)
STANDARD = MachineType(
    'standard',
    ncpu=16,  # 32
    mem_gb_per_core=3.75,  # Total 60G
    price_per_hour=1.0787,  # 32-core machine would cost 2.1574
)

# ~1.7 times memory per core than standard-32, but only 16 cores.
# Total 104G - smaller than standard-32, but less expensive, so useful for 
# memory consuming tools that don't benefit from multiple threads 
HIGHMEM = MachineType(
    'highmem', 
    ncpu=16,
    mem_gb_per_core=6.5,  
    price_per_hour=1.3431,
)


@dataclass(init=False)
class JobResource:
    """
    Represents a fraction of a Hail Batch instance.
    """
    def __init__(
        self,
        machine_type: MachineType,
        ncpu: Optional[int] = None,
        attach_disk_storage_gb: Optional[float] = None,
    ):
        """
        :param machine_type: Hail Batch machine pool type
        :param ncpu: number of CPU request. Will be used to calculate the fraction of
            the machine to take. If not set, all machine's CPUs will be used.
        :param attach_disk_storage_gb: if set to > MachineType.max_default_storage_gb, 
            a larger disc will be attached by Hail Batch. 
        """
        self.machine_type = machine_type
        
        self.fraction_of_full: float = 1.0
        if ncpu is not None:
            if ncpu > self.machine_type.max_ncpu:
                raise ValueError(
                    f'Max number of CPU on machine {self.machine_type.name} '
                    f'is {self.machine_type.max_ncpu}, requested {ncpu}'
                )
            self.fraction_of_full = ncpu / self.machine_type.max_ncpu

        self.attach_disk_storage_gb = None
        if attach_disk_storage_gb is not None:
            if self.fraction_of_full < 1:
                raise ValueError(
                    f'Storage can be overridden only when the entire machine is used, '
                    f'not a fraction ({self.fraction_of_full}). '
                    f'override_storage_gb={attach_disk_storage_gb}'
                )
            self.attach_disk_storage_gb = attach_disk_storage_gb

    def get_mem_gb(self) -> float:
        return self.get_ncpu() * self.machine_type.mem_gb_per_core

    def get_java_mem_mb(self) -> int:
        """
        Calculate memory to pass to the `java -Xms` option. 
        Subtracts 1G to start a java VM, and converts to MB as the option doesn't
        support fractions of GB.
        """
        return int(math.floor(self.get_mem_gb() - 1))

    def get_ncpu(self) -> int:
        return int(self.machine_type.max_ncpu * self.fraction_of_full)
    
    def get_nthreads(self) -> int:
        return self.get_ncpu() * MachineType.THREADS_ON_CPU

    def get_storage_gb(self) -> float:
        """
        Calculate storage in GB
        """
        if self.attach_disk_storage_gb:
            storage_gb = self.attach_disk_storage_gb
        else:
            storage_gb = self.machine_type.calc_instance_disk_gb() * self.fraction_of_full

        # Hail Batch actually requests 5% lower number than the 
        # requested one (e.g. "req_storage: 46.25G, actual_storage: 44.0 GiB"),
        # so we will ask for a bigger number.        
        return storage_gb * 1.05

    def set_to_job(self, j: Job) -> 'JobResource':
        """
        Set the resources to a Job object. Return self to allow chaining, e.g.:
        >>> nthreads = STANDARD.request_resources(nthreads=4).set_to_job(j).get_nthreads()
        """
        j.cpu(self.get_ncpu())
        j.memory(f'{self.get_mem_gb()}G')
        j.storage(f'{self.get_storage_gb()}G')
        # returning self to allow command chaining.
        return self
