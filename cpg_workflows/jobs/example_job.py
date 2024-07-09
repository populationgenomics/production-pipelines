from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import command


def example_job(batch, sequencing_group):
    """
    A test job
    """

    j = batch.new_job('test_job')

    cmd = ''

    if reference_version := config_retrieve('reference_version'):
        # Use a different reference version for your file
        cmd += f'echo "Reference version: {reference_version}"\n'

    if resource_overrides := config_retrieve('resource_override'):
        # Use a different resources for your tools
        cmd += f'echo "Resource override: {resource_overrides}"\n'

    if not config_retrieve('skip_stage', False):
        # Do something
        cmd += f'echo "Hello from {sequencing_group.id}"'

    j.command(command(cmd))

    return j


def example_job(batch, sequencing_group, config):
    """
    A test job
    """

    j = batch.new_job('test_job')

    cmd = ''

    if reference_version := config.reference_version:
        # Use a different reference version for your file
        cmd += f'echo "Reference version: {reference_version}"\n'

    if resource_overrides := config.resource_override:
        # Use a different resources for your tools
        cmd += f'echo "Resource override: {resource_overrides}"\n'

    if not config.skip_stage:
        # Do something
        cmd += f'echo "Hello from {sequencing_group.id}"'

    j.command(command(cmd))

    return j
