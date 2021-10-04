from textwrap import dedent
from typing import Optional, Union, List


def wrap_command(
    command: str,
    output_bucket_path_to_check: Optional[Union[str, List[str]]] = None,
    overwrite: Optional[bool] = False,
    monitor_space: bool = False,
):
    """
    Wraps a command for submission
    If job_resource is defined, monitors output space.
    If output_bucket_path_to_check is defined, checks if this file(s) exists,
    and if it does, skips running the rest of the job.
    """
    return dedent(
        f"""
    set -o pipefail
    set -ex

    {check_existence_command(output_bucket_path_to_check, overwrite)}

    {f'(while true; do {monitor_space_command()}; sleep 600; done) &'
    if monitor_space else ''}

    {command}

    {monitor_space_command() if monitor_space else ''}
    """
    )


def check_existence_command(
    output_path: Optional[Union[str, List[str]]] = None,
    overwrite: bool = True,
) -> str:
    """
    Command that checks the `output_path` existence and exists with rc=0 if it does
    """
    if output_path and not overwrite:
        if isinstance(output_path, str):
            output_path = [output_path]
        return dedent(
            f"""
        # If the output file exists, not running the job
        export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
        gcloud -q auth activate-service-account \\
        --key-file=$GOOGLE_APPLICATION_CREDENTIALS
        gsutil ls {' '.join(output_path)} && exit 0 || echo "Running command"
        """
        )
    return ''


def monitor_space_command():
    """
    Make command that monitors the instance storage space and memory
    """
    return f'df -h; du -sh /io; du -sh /io/batch'
