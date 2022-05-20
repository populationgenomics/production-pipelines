"""
CPG implementations of providers.
"""
import json
import os
from urllib import request

import yaml
from cpg_utils.cloud import read_secret, email_from_id_token

from ... import Namespace
from ...targets import parse_stack

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
    dataset: str | None = None,
    access_level: str | None = None,
):
    """
    Simulate the analysis-runner environment
    """
    access_level = os.getenv('CPG_ACCESS_LEVEL', access_level)
    dataset = os.getenv('CPG_DATASET', dataset)
    if not dataset:
        raise ValueError('CPG_DATASET environment variable must be set')

    namespace = Namespace.TEST if access_level == 'test' else Namespace.MAIN
    stack, namespace = parse_stack(dataset, namespace)
    if namespace == Namespace.TEST or not access_level:
        access_level = 'test'

    sa_email = get_stack_sa_email(stack, access_level)
    if sa_key_template := os.getenv('CPG_SA_PATH_TEMPLATE'):
        sa_key_path = sa_key_template.format(sa_name=sa_email.split('@')[0])
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = sa_key_path

    # overriding in case if dataset was in form of "dataset-test"
    os.environ['CPG_DATASET'] = stack
    os.environ['CPG_ACCESS_LEVEL'] = access_level
    
    hail_token = os.getenv('HAIL_TOKEN')
    gcp_project = os.getenv('CPG_DATASET_GCP_PROJECT')
    if not hail_token or not gcp_project:
        server_config = get_server_config()
        hail_token = server_config[stack][f'{access_level}Token']
        gcp_project = server_config[stack]['projectId']
    assert hail_token
    assert gcp_project

    os.environ.setdefault('CPG_DRIVER_IMAGE', DRIVER_IMAGE)
    os.environ.setdefault('CPG_IMAGE_REGISTRY_PREFIX', IMAGE_REGISTRY_PREFIX)
    os.environ.setdefault('CPG_REFERENCE_PREFIX', REFERENCE_PREFIX)
    os.environ.setdefault('CPG_WEB_URL_TEMPLATE', WEB_URL_TEMPLATE)
    os.environ.setdefault('CPG_OUTPUT_PREFIX', 'cpg-pipes')
    os.environ.setdefault('CPG_DATASET_GCP_PROJECT', gcp_project)
    os.environ.setdefault('HAIL_BILLING_PROJECT', stack)
    os.environ.setdefault('HAIL_BUCKET', f'cpg-{stack}-hail')
    os.environ.setdefault('HAIL_TOKEN', hail_token)
