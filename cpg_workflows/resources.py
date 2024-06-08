"""
Functions to set up Hail Batch resources (cores, memory, storage).
"""

import math
from dataclasses import dataclass

from hailtop.batch.job import Job

from cpg_utils.config import get_config


def _is_power_of_two(n: int) -> bool:
    return math.ceil(math.log(n, 2)) == math.floor(math.log(n, 2))


def gcp_machine_name(name: str, ncpu: int) -> str:
    """
    Machine type name in the GCP world
    """
    assert name in ['standard', 'highmem', 'highcpu'], name
    assert _is_power_of_two(ncpu), ncpu
    return f'n1-{name}-{ncpu}'


@dataclass(init=False)
class MachineType:
    """
    Hail Batch machine type on GCP
    """

    min_cpu: int = 2
    # GCP supports one active thread per vCPU.
    # See https://cloud.google.com/compute/docs/instances/set-threads-per-core
    threads_on_cpu = 1

    def __init__(
        self,
        name: str,
        ncpu: int,
        mem_gb_per_core: float,
        price_per_hour: float,
        disk_size_gb: int,
    ):
        self.name = name
        self.max_ncpu = ncpu
        self.mem_gb_per_core = mem_gb_per_core
        self.price_per_hour = price_per_hour
        self.disk_size_gb = disk_size_gb

    def gcp_name(self):
        """
        Machine type name in the GCP world
        """
        return gcp_machine_name(self.name, self.max_ncpu)

    def max_threads(self) -> int:
        """
        Number of available threads
        """
        return self.max_ncpu * self.threads_on_cpu

    def calc_instance_disk_gb(self) -> int:
        """
        The maximum available storage on an instance is calculated
        in `batch/batch/utils.py/unreserved_worker_data_disk_size_gib()`
        as the disk size (375G) minus reserved image size (30G) minus
        reserved storage per core (5G*ncpu = 120G for a 32-core instance),
        """
        reserved_gb = 30
        reserved_gb_per_core = 5
        return self.disk_size_gb - reserved_gb - reserved_gb_per_core * self.max_ncpu

    def set_resources(
        self,
        j: Job,
        fraction: float | None = None,
        ncpu: int | None = None,
        nthreads: int | None = None,
        mem_gb: float | None = None,
        storage_gb: float | None = None,
    ) -> 'JobResource':
        """
        Set resources to a Job object. If any optional parameters are set,
        they will be used as a bound to request a fraction of an instance.
        """
        return self.request_resources(
            fraction=fraction,
            ncpu=ncpu,
            nthreads=nthreads,
            mem_gb=mem_gb,
            storage_gb=storage_gb,
        ).set_to_job(j)

    def request_resources(
        self,
        fraction: float | None = None,
        ncpu: int | None = None,
        nthreads: int | None = None,
        mem_gb: float | None = None,
        storage_gb: float | None = None,
    ) -> 'JobResource':
        """
        Request resources from the machine, satisfying all provided requirements.
        If not requirements are provided, the minimal amount of cores
        (self.MIN_NCPU) will be used.
        """
        # determining the biggest limit to satisfy, measured in the number of CPUs:
        min_ncpu = max(
            filter(
                None,
                [
                    self.adjust_ncpu(ncpu or self.min_cpu),
                    self.fraction_to_ncpu(fraction) if fraction else None,
                    self.nthreads_to_ncpu(nthreads) if nthreads else None,
                    self.mem_gb_to_ncpu(mem_gb) if mem_gb else None,
                    self.storage_gb_to_ncpu(storage_gb) if storage_gb else None,
                ],
            ),
        )
        return JobResource(
            machine_type=self,
            ncpu=min_ncpu,
            attach_disk_storage_gb=(storage_gb if storage_gb and storage_gb > self.calc_instance_disk_gb() else None),
        )

    def fraction_to_ncpu(self, fraction: float) -> int:
        """
        Converts fraction to the number of CPU (e.g. fraction=1.0 to take the entire
        machine, fraction=0.5 to take half of it, etc.).
        """
        ncpu = int(math.ceil(self.max_ncpu * fraction))
        return self.adjust_ncpu(ncpu)

    def mem_gb_to_ncpu(self, mem_gb: float) -> int:
        """
        Converts memory requirement to the number of CPU requirement.
        """
        ncpu = int(math.ceil(mem_gb / self.mem_gb_per_core))
        return self.adjust_ncpu(ncpu)

    def storage_gb_to_ncpu(self, storage_gb: float) -> int:
        """
        Converts storage requirement to the number of CPU requirement.

        We want to avoid attaching disks: attaching a disk to an existing instance
        might fail with `mkfs.ext4 ...` error, see:
        https://batch.hail.populationgenomics.org.au/batches/7488/jobs/12
        So this function will calculate the number of CPU to request so your jobs
        can be packed to fit the default instance's available storage
        (calculated with self.calc_instance_disk_gb()).
        """
        fraction = storage_gb / self.calc_instance_disk_gb()
        fraction = min(fraction, 1.0)
        return self.fraction_to_ncpu(fraction)

    def nthreads_to_ncpu(self, nthreads: int) -> int:
        """
        Convert number of threads into number of cores/CPU
        """
        return self.adjust_ncpu(math.ceil(nthreads / 2))

    def adjust_ncpu(self, ncpu: int) -> int:
        """
        Adjust request number of CPU to a number allowed by Hail, i.e.
        the nearest power of 2, not less than the minimal number of cores allowed.
        """
        if ncpu > self.max_ncpu:
            raise ValueError(f'Requesting more cores than available on {self.name} machine: {ncpu}>{self.max_ncpu}')

        if ncpu < MachineType.min_cpu:
            ncpu = MachineType.min_cpu

        # round to the nearest power of 2 (15 -> 16, 16 -> 16, 17 -> 32)
        return int(pow(2, math.ceil(math.log2(ncpu))))


# Default Hail Batch machine.
#
# Bigger default number of cores (32 vs 16) would allow for more threads for alignment,
# however would mean a smaller available storage (5G per core is reserved)
STANDARD = MachineType(
    'standard',
    ncpu=16,
    mem_gb_per_core=3.75,  # Total 60G
    price_per_hour=1.0787,  # 32-core machine would cost 2.1574
    disk_size_gb=375,
)

# ~1.7 times memory per core than standard-32, but only 16 cores.
# Total 104G - smaller than standard-32, but less expensive, so useful for
# memory consuming tools that don't benefit from multiple threads
HIGHMEM = MachineType(
    'highmem',
    ncpu=16,
    mem_gb_per_core=6.5,
    price_per_hour=1.3431,
    disk_size_gb=375,
)


@dataclass(init=False)
class JobResource:
    """
    Represents a fraction of a Hail Batch instance.
    """

    def __init__(
        self,
        machine_type: MachineType,
        ncpu: int | None = None,
        attach_disk_storage_gb: float | None = None,
    ):
        """
        @param machine_type: Hail Batch machine pool type
        @param ncpu: number of CPU request. Will be used to calculate the fraction of
            the machine to take. If not set, all machine's CPUs will be used.
        @param attach_disk_storage_gb: if set to > MachineType.max_default_storage_gb,
            a larger disc will be attached by Hail Batch.
        """
        self.machine_type = machine_type

        self.fraction_of_full: float = 1.0
        if ncpu is not None:
            if ncpu > self.machine_type.max_ncpu:
                raise ValueError(
                    f'Max number of CPU on machine {self.machine_type.name} '
                    f'is {self.machine_type.max_ncpu}, requested {ncpu}',
                )
            self.fraction_of_full = ncpu / self.machine_type.max_ncpu

        self.attach_disk_storage_gb = None
        if attach_disk_storage_gb is not None:
            if self.fraction_of_full < 1:
                raise ValueError(
                    f'Storage can be overridden only when the entire machine is used, '
                    f'not a fraction ({self.fraction_of_full}). '
                    f'override_storage_gb={attach_disk_storage_gb}',
                )
            self.attach_disk_storage_gb = attach_disk_storage_gb

    def get_mem_gb(self) -> float:
        """
        Memory resources in GB
        """
        return self.get_ncpu() * self.machine_type.mem_gb_per_core

    def java_mem_options(self, overhead_gb: float = 1) -> str:
        """
        Returns -Xms -Xmx options to set Java JVM memory usage to use all the memory
        resources represented.
        @param overhead_gb: Amount of memory (in decimal GB) to leave available for
        other purposes.
        """
        mem_bytes = (self.get_mem_gb() - overhead_gb) * 1_000_000_000
        # Approximate as binary MiB (but not GiB as these options don't support
        # fractional values) so that logs are easier to read
        mem_mib = math.floor(mem_bytes / 1_048_576)
        return f'-Xms{mem_mib}M -Xmx{mem_mib}M'

    def java_gc_thread_options(self, surplus: int = 2) -> str:
        """
        Returns -XX options to set Java JVM garbage collection threading.
        @param surplus: Number of threads to leave available for other purposes.
        """
        gc_threads = self.get_nthreads() - surplus
        return f'-XX:+UseParallelGC -XX:ParallelGCThreads={gc_threads}'

    def get_ncpu(self) -> int:
        """
        Number of cores/CPU
        """
        return int(self.machine_type.max_ncpu * self.fraction_of_full)

    def get_nthreads(self) -> int:
        """
        Number of threads
        """
        return self.get_ncpu() * MachineType.threads_on_cpu

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

        j.storage(f'{self.get_storage_gb()}G')
        j.cpu(self.get_ncpu())
        j.memory(f'{self.get_mem_gb()}G')

        # Returning self to allow command chaining.
        return self


def storage_for_cram_qc_job() -> int | None:
    """
    Get storage request for a CRAM QC processing job, gb
    """
    sequencing_type = get_config()['workflow']['sequencing_type']
    storage_gb = None  # avoid extra disk by default
    if sequencing_type == 'genome':
        storage_gb = 100
    if sequencing_type == 'exome':
        storage_gb = 20
    return storage_gb


def joint_calling_scatter_count(sequencing_group_count: int) -> int:
    """
    Number of partitions for joint-calling jobs (GenotypeGVCFs, VQSR, VEP),
    as a function of the sequencing group number.
    """
    if scatter_count := get_config()['workflow'].get('scatter_count'):
        return scatter_count

    # Estimating this is challenging because GenotypeGVCFs does not scale
    # linearly with the number of genomes.
    # Values are adjusted based on experience with the actual number of genomes.
    # e.g. 1000 scatter count was too low for 3800 genomes.
    for threshold, scatter_count in {
        4000: 1400,
        3500: 1200,
        3000: 1000,
        2000: 600,
        1000: 400,
        500: 200,
        250: 100,
    }.items():
        if sequencing_group_count >= threshold:
            return scatter_count
    return 50


def storage_for_joint_vcf(
    sequencing_group_count: int | None,
    site_only: bool = True,
) -> float | None:
    """
    Storage enough to fit and process a joint-called VCF
    """
    if not sequencing_group_count:
        return None
    if get_config()['workflow']['sequencing_type'] == 'exome':
        gb_per_sequencing_group = 0.1
    else:
        gb_per_sequencing_group = 1.0
        if not site_only:
            gb_per_sequencing_group = 1.5

    return gb_per_sequencing_group * sequencing_group_count
