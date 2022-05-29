"""
CPG implementations of providers.
"""
import json
import os
import tempfile
import uuid
from urllib import request

import toml  # noqa
import yaml
from cpg_utils.cloud import read_secret
from cpg_utils.config import set_config_path

from ... import Namespace, to_path, Path
from ...targets import parse_stack
from ...utils import exists

ANALYSIS_RUNNER_PROJECT_ID = 'analysis-runner'

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:7f41eeb90c8bec8836a1cd20ad1911b5989a5893-hail-db65c33c29100c64405c39ebace90a7c463b4bec'
IMAGE_REGISTRY_PREFIX = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
REFERENCE_PREFIX = 'gs://cpg-reference'
WEB_URL_TEMPLATE = f'https://{{namespace}}-web.populationgenomics.org.au/{{dataset}}'


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


def analysis_runner_environment(
    local_tmp_dir: Path | None = None,
    dataset: str | None = None,
    namespace: Namespace | None = None,
):
    """
    Simulate the analysis-runner environment. Requires only the CPG_DATASET environment 
    variable or `dataset` parameter set. Optionally, SERVICE_ACCOUNT_KEY can be set
    to run as a specific user. In this case, GOOGLE_APPLICATION_CREDENTIALS should be
    set with permissions to pull Hail tokens secret.
    
    Writes a config to gs://cpg-config and sets CPG_CONFIG_PATH, the only required
    environment variable for cpg-utils. Also sets HAIL_TOKEN.
    """
    if not local_tmp_dir:
        local_tmp_dir = to_path(tempfile.mkdtemp())
    assert local_tmp_dir

    dataset = os.getenv('CPG_DATASET', dataset)
    if not dataset:
        raise ValueError(
            'CPG_DATASET environment variable or dataset parameter must be set'
        )

    stack, namespace = parse_stack(dataset, namespace)
    access_level = 'test' if namespace == Namespace.TEST else 'full'

    # Parse server config to get Hail token and the GCP project name
    hail_token = os.getenv('HAIL_TOKEN')
    gcp_project = os.getenv('CPG_DATASET_GCP_PROJECT')
    if not hail_token or not gcp_project:
        server_config = get_server_config()
        hail_token = server_config[stack][f'{access_level}Token']
        gcp_project = server_config[stack]['projectId']
    assert hail_token
    assert gcp_project
    os.environ['HAIL_TOKEN'] = hail_token
    
    # Active dataset service account
    if sa_key_path := os.getenv('SERVICE_ACCOUNT_KEY'):
        if '{name}' in sa_key_path:
            sa_email = get_stack_sa_email(stack, access_level)
            sa_key_path = sa_key_path.format(name=sa_email.split('@')[0])
        if not exists(sa_key_path):
            raise ValueError(f'Service account key does not exist: {sa_key_path}')
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = sa_key_path

    config = {
        'hail': {
            'billing_project': dataset,
            'bucket': f'cpg-{stack}-hail',
        },
        'workflow': {
            'access_level': access_level,
            'dataset': stack,
            'dataset_gcp_project': gcp_project,
            'driver_image': DRIVER_IMAGE,
            'image_registry_prefix': IMAGE_REGISTRY_PREFIX,
            'reference_prefix': REFERENCE_PREFIX,
            'output_prefix': 'cpg-pipes',
            'web_url_template': WEB_URL_TEMPLATE,
        },
    }

    config_path = local_tmp_dir / (str(uuid.uuid4()) + '.toml')
    with config_path.open('w') as f:
        toml.dump(config, f)

    set_config_path(config_path)
