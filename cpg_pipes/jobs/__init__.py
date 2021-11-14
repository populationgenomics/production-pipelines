def wrap_command(
    command: str,
    monitor_space: bool = False,
    setup_gcp: bool = False,
    dedent: bool = True
) -> str:
    """
    Wraps a command for submission
    If job_resource is defined, monitors output space.
    If output_bucket_path_to_check is defined, checks if this file(s) exists,
    and if it does, skips running the rest of the job.
    """
    gcp_cmd = ''
    if setup_gcp:
        gcp_cmd = """\
        export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
        gcloud -q auth activate-service-account \
        --key-file=$GOOGLE_APPLICATION_CREDENTIALS
        """
    
    cmd = f"""\
    set -o pipefail
    set -ex
    {gcp_cmd}
    
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


# def check_existence_command(
#     output_path: Optional[Union[str, List[str]]] = None,
#     overwrite: bool = True,
# ) -> str:
#     """
#     Command that checks the `output_path` existence and exists with rc=0 if it does
#     """
#     if output_path and not overwrite:
#         if isinstance(output_path, str):
#             output_path = [output_path]
#         return dedent(f"""\
#         # If the output file exists, not running the job
#         export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
#         gcloud -q auth activate-service-account \\
#         --key-file=$GOOGLE_APPLICATION_CREDENTIALS
#         gsutil ls {' '.join(output_path)} && (exit 0) || echo "Running command"
#         """)
#     return ''


def monitor_space_command():
    """
    Make command that monitors the instance storage space and memory
    """
    return f'df -h; du -sh /io; du -sh /io/batch'


# def new_job(
#     b: hb.Batch, 
#     name: str, 
#     sample_name: Optional[str] = None, 
#     project_name: Optional[str] = None,
#     **kwargs
# ) -> Job:
#     """
#     Simplifies calling a new job by wrapping  arguments into a
#     attributes dictionary
#     """
#     assert isinstance(b, hb.Batch)
#     if sample_name:
#         kwargs['sample'] = sample_name
#     if project_name:
#         kwargs['project'] = project_name
#     return b.new_job(name, attributes=kwargs)
