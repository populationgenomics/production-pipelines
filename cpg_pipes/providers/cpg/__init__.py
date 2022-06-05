"""
CPG implementations of providers.
"""
import json
import os
import tempfile
import uuid
from urllib import request

import toml
import yaml
from cpg_utils.cloud import read_secret
from cpg_utils.config import set_config_path, update_dict, get_config

from ... import Namespace, to_path, Path
from ...targets import parse_stack
from ...utils import exists

ANALYSIS_RUNNER_PROJECT_ID = 'analysis-runner'


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


def overwrite_config(config: dict) -> dict:
    """
    Writes and sets cpg-utils config.
    """
    get_config()['workflow'].setdefault('local_tmp_dir', tempfile.mkdtemp())
    local_tmp_dir = to_path(get_config()['workflow']['local_tmp_dir'])
    config_path = local_tmp_dir / (str(uuid.uuid4()) + '.toml')
    with config_path.open('w') as f:
        toml.dump(config, f)
    set_config_path(config_path)
    return get_config()


def complete_infra_config(config: dict) -> dict:
    """
    Set up the infrastructure config parameters, unless they are already set
    by analysis-runner.

    Note for CPG: if typical analysis-runner configuration is not included, this 
    function would somulate the analysis-runner environment. 
    It requires only the `dataset` parameter, however optionally SERVICE_ACCOUNT_KEY 
    can be set to impersonate specific service account. In this case, 
    GOOGLE_APPLICATION_CREDENTIALS should also initially be set to the analysis-runner 
    service account to allow read permissions to the secret with Hail tokens.
    """
    dataset = config['workflow']['dataset']
    access_level = config['workflow'].get('access_level', 'standard')
    dataset_gcp_project = config['workflow'].get('dataset_gcp_project')
    
    stack, namespace = parse_stack(dataset, Namespace.from_access_level(access_level))
    if not access_level:
        access_level = 'test' if namespace == Namespace.TEST else 'full'
        config['workflow'].setdefault('access_level', access_level)

    # Parse server config to get Hail token and the GCP project name
    hail_token = os.getenv('HAIL_TOKEN')
    if not hail_token or not dataset_gcp_project:
        server_config = get_server_config()
        hail_token = server_config[stack][f'{access_level}Token']
        dataset_gcp_project = server_config[stack]['projectId']
        config['workflow'].setdefault('dataset_gcp_project', dataset_gcp_project)

    assert hail_token
    assert dataset_gcp_project
    os.environ['HAIL_TOKEN'] = hail_token
    
    # Active dataset service account
    if sa_key_path := os.getenv('SERVICE_ACCOUNT_KEY'):
        if '{name}' in sa_key_path:
            sa_email = get_stack_sa_email(stack, access_level)
            sa_key_path = sa_key_path.format(name=sa_email.split('@')[0])
        if not exists(sa_key_path):
            raise ValueError(f'Service account key does not exist: {sa_key_path}')
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = sa_key_path
    
    config['workflow'].setdefault('dataset', stack)
    config['workflow'].setdefault('output_prefix', 'cpg-pipes')
    config['hail'].setdefault('billing_project', stack)
    config['hail'].setdefault('bucket', f'cpg-{stack}-hail')
    return config
