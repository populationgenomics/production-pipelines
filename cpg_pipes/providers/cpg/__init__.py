"""
CPG implementations of providers.
"""
import json
import os

from cpg_utils.cloud import read_secret

from ... import Namespace
from ...targets import parse_stack

ANALYSIS_RUNNER_PROJECT_ID = 'analysis-runner'

DRIVER_IMAGE = 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:7f41eeb90c8bec8836a1cd20ad1911b5989a5893-hail-db65c33c29100c64405c39ebace90a7c463b4bec'
IMAGE_REGISTRY_PREFIX = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
REFERENCE_PREFIX = 'gs://cpg-reference'


def get_server_config() -> dict:
    """Get the server-config from the secret manager"""
    return json.loads(read_secret(ANALYSIS_RUNNER_PROJECT_ID, 'server-config'))


def analysis_runner_environment(
    dataset: str | None = None,
    namespace: Namespace | None = None,
):
    """Simulate analysis-runner environment"""
    dataset = dataset or os.environ['ANALYSIS_DATASET']
    stack, namespace = parse_stack(dataset, namespace)

    server_config = get_server_config()
    
    access_level = 'test' if namespace == Namespace.TEST else 'full'
    hail_token = server_config[stack].get(f'{access_level}Token')

    os.environ['CPG_ACCESS_LEVEL'] = access_level
    os.environ['CPG_DATASET'] = stack
    os.environ['CPG_DATASET_GCP_PROJECT'] = server_config[stack]['projectId']
    os.environ['CPG_DRIVER_IMAGE'] = DRIVER_IMAGE
    os.environ['CPG_IMAGE_REGISTRY_PREFIX'] = IMAGE_REGISTRY_PREFIX
    os.environ['CPG_REFERENCE_PREFIX'] = REFERENCE_PREFIX
    os.environ['CPG_OUTPUT_PREFIX'] = 'cpg-pipes'
    os.environ['HAIL_BILLING_PROJECT'] = stack
    os.environ['HAIL_BUCKET'] = f'cpg-{stack}-hail'
    os.environ['HAIL_TOKEN'] = hail_token


analysis_runner_environment()
