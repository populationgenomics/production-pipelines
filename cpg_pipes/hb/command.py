"""
Helpers to setup Job's command.
"""

import logging
from typing import List, Union


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

# command that monitors the instance storage space
MONITOR_SPACE_CMD = f'df -h; du -sh /io; du -sh /io/batch'


def wrap_command(
    command: Union[str, List[str]],
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

    @param command: command to wrap (can be a list of commands)
    @param monitor_space: add a background process that checks the instance disk 
        space every 5 minutes and prints it to the screen
    @param setup_gcp: login to GCP
    @param define_retry_function: when set, adds bash functions `retry` that attempts 
        to redo a command with a pause of default 30 seconds (useful to pull inputs 
        and get around GoogleEgressBandwidth Quota or other google quotas)
    @param dedent: remove all common leading intendation from the command
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
    
    {command}
    
    {MONITOR_SPACE_CMD if monitor_space else ''}
    """

    if dedent:
        # remove any leading spaces and tabs
        cmd = '\n'.join(line.strip() for line in cmd.split('\n'))
        # remove sretches of spaces
        cmd = '\n'.join(' '.join(line.split()) for line in cmd.split('\n'))
    return cmd
