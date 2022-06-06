"""
CPG implementations of providers.
"""
import json
import logging
import os
import tempfile
import uuid
from contextlib import contextmanager
from urllib import request

import toml
import yaml
from cloudpathlib import AnyPath
from cpg_utils.cloud import read_secret
from cpg_utils.config import set_config_paths, get_config

from ...utils import exists

ANALYSIS_RUNNER_PROJECT_ID = 'analysis-runner'

logger = logging.getLogger(__file__)


def get_server_config() -> dict:
    """Get the server-config from the secret manager"""
    secret = read_secret(ANALYSIS_RUNNER_PROJECT_ID, 'server-config')
    if not secret:
        raise ValueError(
            f'Cannot read secret "server-config" from the GCP project '
            f'"{ANALYSIS_RUNNER_PROJECT_ID}". Please, set HAIL_TOKEN and '
            f'CPG_DATASET_GCP_PROJECT environment variables explicitly.'
        )
    return json.loads(secret)


def get_stack_sa_email(stack: str, access_level: str) -> str:
    """Parse service account email from Pulumi yaml"""
    seqr_stack_url = (
        f'https://raw.githubusercontent.com/populationgenomics/analysis-runner/main'
        f'/stack/Pulumi.{stack}.yaml'
    )
    with request.urlopen(seqr_stack_url) as f:
        data = yaml.safe_load(f)['config']
    return data[f'datasets:hail_service_account_{access_level}']


@contextmanager
def analysis_runner_env():
    """
    Context manager version of `set_analysis_runner_env`. Will revert the config back 
    after existing.
    """
    # Keeping the original environment to revert later:
    original_config_path = os.getenv('CPG_CONFIG_PATH')
    original_sa_key_path = os.getenv('GOOGLE_APPLICATION_CREDENTIALS')
    try:
        if sa_key_path := os.getenv('SERVICE_ACCOUNT_KEY'):
            set_analysis_runner_env(sa_key_path)
        yield
    finally:  # Revert the environment back to the original
        set_config_paths([original_config_path])
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = original_sa_key_path


def set_analysis_runner_env(sa_key_path: str):
    """
    @param sa_key_path: path to the dataset service account. Can be parametrised
    by service account name in curly braces, e.g. "/path/to/{name}_key.json".

    Set up the infrastructure config parameters. Requires 
    GOOGLE_APPLICATION_CREDENTIALS to be set to the analysis-runner service account 
    to allow read permissions to the secret with Hail tokens (unless HAIL_TOKEN and 
    `workflow.dataset_gcp_project` are already set).

    Required config fields:
     - `workflow.dataset`
     - `workflow.access_level`

    Config parameters to be added:
    - `workflow.output_prefix`
    - `workflow.dataset_gcp_project`
    - `workflow.access_level`
    - `hail.billing_project`
    - `hail.bucket`

    Also sets the `HAIL_TOKEN` environment variable.
    """
    if not os.getenv('GOOGLE_APPLICATION_CREDENTIALS'):
        raise ValueError(
            'GOOGLE_APPLICATION_CREDENTIALS must be set and point to the '
            'analysis-runner service account'
        )

    dataset = get_config()['workflow']['dataset']
    access_level = get_config()['workflow']['access_level']

    # Parse server config to get Hail token and the GCP project name
    hail_token = os.getenv('HAIL_TOKEN')
    dataset_gcp_project = get_config()['workflow'].get('dataset_gcp_project')
    if not hail_token or not dataset_gcp_project:
        server_config = get_server_config()
        hail_token = server_config[dataset][f'{access_level}Token']
        dataset_gcp_project = server_config[dataset]['projectId']
        assert hail_token  # assert to make type checking happy
        get_config()['workflow'].setdefault('dataset_gcp_project', dataset_gcp_project)
        os.environ.setdefault('HAIL_TOKEN', hail_token)

    # Activate dataset service account
    if '{name}' in sa_key_path:  # Path can be parametrised by a SA name.
        sa_email = get_stack_sa_email(dataset, access_level)
        sa_key_path = sa_key_path.format(name=sa_email.split('@')[0])
    if not exists(sa_key_path):
        raise ValueError(f'Service account key does not exist: {sa_key_path}')
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = sa_key_path

    get_config()['workflow'].setdefault('output_prefix', 'cpg-pipes')
    get_config()['hail'].setdefault('billing_project', dataset)
    get_config()['hail'].setdefault('bucket', f'cpg-{dataset}-hail')
    tmp_dir = get_config()['workflow'].setdefault('local_tmp_dir', tempfile.mkdtemp())
    write_config(get_config(), tmp_dir)


def write_config(config: dict, prefix: str, path: str | None = None) -> None:
    """Writes config to a file, and sets the path to that file as a new config path.

    Parameters
    ----------
    config: dict
        Config object
    prefix: str
        Path prefix (directory) where the TOML file will be written
    path: str, optional
        Path to write the config, takes precedence over prefix
    """
    if not path:
        path = os.path.join(prefix, f'{uuid.uuid4()}.toml')
    with AnyPath(path).open('w') as f:
        toml.dump(config, f)
    set_config_paths([path])
