"""
Helpers to set up Job's command.
"""
import inspect
import logging
from typing import List, Union
import textwrap

from cpg_utils.config import get_config
from cpg_utils.hail_batch import genome_build, remote_tmpdir
from hailtop.batch import ResourceFile

from cpg_pipes import Path

logger = logging.getLogger(__file__)


# commands that activate gsutil
GCLOUD_CMD = """\
export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
gcloud -q auth activate-service-account \
--key-file=$GOOGLE_APPLICATION_CREDENTIALS
"""

# commands that declare functions that pull files on an instance,
# handling transitive errors
RETRY_CMD = """\
function fail {
  echo $1 >&2
  exit 1
}

function retry {
  local n_attempts=10
  local delay=30
  local n=1
  while ! eval "$@"; do
    if [[ $n -lt $n_attempts ]]; then
      ((n++))
      echo "Command failed. Attempt $n/$n_attempts after ${delay}s..."
      sleep $delay;
    else
      fail "The command has failed after $n attempts."
    fi
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

# command that monitors the instance storage space
MONITOR_SPACE_CMD = f'df -h; du -sh /io; du -sh /io/batch'

ADD_SCRIPT_CMD = """\
cat <<EOT >> {script_name}
{script_contents}
EOT\
"""


def python_command(
    module,
    func_name: str,
    *func_args,
    setup_gcp: bool = False,
    setup_hail: bool = False,
    packages: list[str] | None = None,
):
    """
    Construct a command for a Job that runs a python function.
    If hail_billing_project is provided, Hail Query will be also initialised.
    """
    billing_project = get_config()['hail']['billing_project']
    dataset = get_config()['workflow']['dataset']
    bucket = remote_tmpdir(f'cpg-{dataset}-hail')

    python_cmd = f"""
import logging
logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)
"""
    if setup_hail:
        python_cmd += f"""
import asyncio
import hail as hl
asyncio.get_event_loop().run_until_complete(
    hl.init_batch(
        default_reference='{genome_build()}',
        billing_project='{billing_project}',
        remote_tmpdir='{bucket}',
    )
)
"""
    python_cmd += f"""
{textwrap.dedent(inspect.getsource(module))}
{func_name}{func_args}
"""
    cmd = f"""
set -o pipefail
set -ex
{GCLOUD_CMD if setup_gcp else ''}

{('pip3 install ' + ' '.join(packages)) if packages else ''}

cat << EOT >> script.py
{python_cmd}
EOT
python3 script.py
"""
    return cmd


def wrap_command(
    command: Union[str, List[str]],
    monitor_space: bool = False,
    setup_gcp: bool = False,
    define_retry_function: bool = False,
    rm_leading_space: bool = True,
    python_script_path: Path | None = None,
) -> str:
    """
    Wraps a command for submission
    If job_resource is defined, monitors output space.
    If output_bucket_path_to_check is defined, checks if this file(s) exists,
    and if it does, skips running the rest of the job.

    @param command: command to wrap (can be a list of commands)
    @param monitor_space: add a background process that checks the instance disk
        space every 5 minutes and prints it to the screen
    @param setup_gcp: login to GCP
    @param define_retry_function: when set, adds bash functions `retry` that attempts
        to redo a command with a pause of default 30 seconds (useful to pull inputs
        and get around GoogleEgressBandwidth Quota or other google quotas)
    @param rm_leading_space: remove all leading spaces and tabs from the command lines
    @param python_script_path: if provided, copy this python script
    """
    if isinstance(command, list):
        command = '\n'.join(command)

    cmd = f"""\
    set -o pipefail
    set -ex
    {GCLOUD_CMD if setup_gcp else ''}
    {RETRY_CMD if define_retry_function else ''}
    
    {f'(while true; do {MONITOR_SPACE_CMD}; sleep 600; done) &'
    if monitor_space else ''}

    {{copy_script_cmd}}

    {command}
    
    {MONITOR_SPACE_CMD if monitor_space else ''}
    """

    if rm_leading_space:
        # remove any leading spaces and tabs
        cmd = '\n'.join(line.strip() for line in cmd.split('\n'))
        # remove sretches of spaces
        cmd = '\n'.join(' '.join(line.split()) for line in cmd.split('\n'))
    else:
        # Remove only common leading space:
        cmd = textwrap.dedent(cmd)

    # We don't want the python script tabs to be stripped, so
    # we are inserting it after leadings space is removed
    if python_script_path:
        with python_script_path.open() as f:
            script_contents = f.read()
        cmd = cmd.replace(
            '{copy_script_cmd}',
            ADD_SCRIPT_CMD.format(
                script_name=python_script_path.name,
                script_contents=script_contents,
            ),
        )
    else:
        cmd = cmd.replace('{copy_script_cmd}', '')

    return cmd


def seds_to_extend_sample_ids(
    rich_id_map: dict[str, str],
    fnames: list[str | ResourceFile],
) -> str:
    """
    Helper function to add seds into a command that would extend samples IDs
    in each file in `fnames` with an external ID, only if external ID is
    different from the original.

    @param rich_id_map: map used to replace samples, e.g. {'CPG1': 'CPG1|MYID'}
    @param fnames: file names and Hail Batch Resource files where to replace IDs
    @return: bash command that does replacement
    """
    cmd = ''
    for sid, rich_sid in rich_id_map.items():
        for fname in fnames:
            cmd += f'sed -iBAK \'s/{sid}/{rich_sid}/g\' {fname}'
            cmd += '\n'
    return cmd
