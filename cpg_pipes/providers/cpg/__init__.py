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
from cpg_utils.config import set_config_path

from ... import Namespace, to_path, Path
from ...targets import parse_stack
from ...utils import exists

ANALYSIS_RUNNER_PROJECT_ID = 'analysis-runner'

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


def set_config(config: dict, local_tmp_dir: Path | None = None):
    """
    Writes and sets cpg-utils config.
    """
    local_tmp_dir = local_tmp_dir or to_path(tempfile.mkdtemp())
    config_path = local_tmp_dir / (str(uuid.uuid4()) + '.toml')
    with config_path.open('w') as f:
        toml.dump(config, f)
    set_config_path(config_path)


def build_infra_config(
    dataset: str,
    namespace: Namespace | None = None,
    gcp_project: str | None = None,
    image_registry_prefix: str = IMAGE_REGISTRY_PREFIX,
    reference_prefix: str = REFERENCE_PREFIX,
    web_url_template: str = WEB_URL_TEMPLATE,
) -> dict[str, dict[str, str]]:
    """
    Set up the infrastructure config.
    
    Note for CPG: if typical analysis-runner configuration is not included, this 
    function would somulate the analysis-runner environment. 
    It requires only the `dataset` parameter, however optionally SERVICE_ACCOUNT_KEY 
    can be set to impersonate specific service account. In this case, 
    GOOGLE_APPLICATION_CREDENTIALS should also initially be set to the analysis-runner 
    service account to allow read permissions to the secret with Hail tokens.
    """
    stack, namespace = parse_stack(dataset, namespace)
    access_level = 'test' if namespace == Namespace.TEST else 'full'

    # Parse server config to get Hail token and the GCP project name
    hail_token = os.getenv('HAIL_TOKEN')
    gcp_project = gcp_project or os.getenv('CPG_DATASET_GCP_PROJECT')
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

    return {
        'hail': {
            'billing_project': dataset,
            'bucket': f'cpg-{stack}-hail',
        },
        'workflow': {
            'access_level': access_level,
            'dataset': stack,
            'dataset_gcp_project': gcp_project,
            'image_registry_prefix': image_registry_prefix or IMAGE_REGISTRY_PREFIX,
            'reference_prefix': reference_prefix or REFERENCE_PREFIX,
            'web_url_template': web_url_template or WEB_URL_TEMPLATE,
            'output_prefix': 'cpg-pipes',
        },
    }
