"""
Perform aberrant splicing analysis with FRASER.
"""

import hailtop.batch as hb
from hailtop.batch.job import Job
from cpg_utils import Path, to_path
from cpg_utils.hail_batch import command, image_path
from cpg_utils.config import get_config
from cpg_workflows.utils import can_reuse
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import (
    BamPath,
)
from cpg_workflows.workflow import (
    SequencingGroup,
)
from textwrap import dedent
from os.path import basename


class Fraser:
    """
    Construct a FRASER command for performing aberrant splicing analysis.
    """

    def __init__(
            self,
            input_bams: list[str | Path],
            txdb_file: str | Path | None = None,
            output_tar_gz_path: str | Path | None = None,
            nthreads: int = 8,
        ) -> None:
        self.input_bams = input_bams
        assert isinstance(self.input_bams, list), f'input_bams must be a list, instead got {self.input_bams}'
        self.input_bams_r_str = ', '.join([f'"{str(f)}"' for f in self.input_bams])
        self.txdb_file_path = str(txdb_file)
        self.output_tar_gz_path = str(output_tar_gz_path)
        self.nthreads = nthreads

        # Build OUTRIDER command
        self.command = """
        R --vanilla <<EOF
        library(FRASER)
        EOF
        """
        # Tar up outputs
        self.command += f"""
        tar -czvf {self.output_tar_gz_path} plots results fraser.results.RData
        """
        self.command = dedent(self.command).strip()

    def __str__(self):
        return self.command
    
    def __repr__(self):
        return self.__str__()


def fraser(
    b: hb.Batch,
    input_bams: list[BamPath | hb.ResourceGroup],
    output_path: str | Path | None = None,
    cohort_name: str | None = None,
    job_attrs: dict[str, str] | None = None,
    requested_nthreads: int | None = None,
) -> Job:
    """
    Run FRASER.
    """

    # Localise input files
    assert all([isinstance(f, (BamPath, hb.ResourceGroup)) for f in input_bams])
    input_bams_localised = [
        f if isinstance(f, hb.ResourceGroup) else b.read_input_group(**{
            'bam': str(f.path),
            'bam.bai': str(f.index_path),
        })
        for f in input_bams
    ]
    assert all([isinstance(hb.ResourceGroup) for f in input_bams_localised])

    # TODO: complete job. Everything below this is copied from outrider job

    gtf_file = get_config()['references'].get('gtf')
    gtf_file = to_path(gtf_file)
    gtf_file_rg = b.read_input_group(gtf=str(gtf_file))

    # Create job
    job_name = f'outrider_{cohort_name}' if cohort_name else 'count'
    _job_attrs = (job_attrs or {}) | dict(label=job_name, tool='outrider')
    j = b.new_job(job_name, _job_attrs)
    # j.image(image_path('outrider'))
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/outrider:1.18.1')

    # Set resource requirements
    nthreads = requested_nthreads or 8
    res = STANDARD.set_resources(
        j,
        ncpu=nthreads,
        storage_gb=50,
    )

    j.declare_resource_group(
        output={
            'tar_gz': '{root}.tar.gz'
        }
    )

    # Create counting command
    outrider = Outrider(
        input_counts=infiles_localised,
        gtf_file=str(gtf_file_rg.gtf),
        output_tar_gz_path=j.output.tar_gz,
        nthreads=res.get_nthreads(),
    )

    cmd = str(outrider)
    j.command(command(cmd, monitor_space=True))

    # Write output to file
    if output_path:
        # NOTE: j.output is just a placeholder
        b.write_output(j.output.tar_gz, str(output_path))
    
    return j
