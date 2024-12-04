import json
from typing import Literal

from google.cloud import secretmanager


def get_ica_secrets() -> dict[Literal['projectID', 'apiKey'], str]:
    SECRET_CLIENT = secretmanager.SecretManagerServiceClient()
    SECRET_PROJECT = 'cpg-common'
    SECRET_NAME = 'illumina_cpg_workbench_api'
    SECRET_VERSION = 'latest'
    try:
        secret_path: str = SECRET_CLIENT.secret_version_path(
            project=SECRET_PROJECT,
            secret=SECRET_NAME,
            secret_version=SECRET_VERSION,
        )
        response: secretmanager.AccessSecretVersionResponse = SECRET_CLIENT.access_secret_version(
            request={'name': secret_path},
        )
        return json.loads(response.payload.data.decode('UTF-8'))
    except Exception as e:
        raise Exception(f'Could not obtain ICA credentials: {e}') from e
