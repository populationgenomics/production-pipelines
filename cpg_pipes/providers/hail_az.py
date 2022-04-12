"""
Extending cloudpathlib Azure implementation to support hail-az:// and https:// schemes.
Inspired by https://github.com/drivendataorg/cloudpathlib/issues/157
"""

import re
from urllib.parse import urlparse

from cloudpathlib import AzureBlobClient, AzureBlobPath
from cloudpathlib.client import register_client_class
from cloudpathlib.cloudpath import register_path_class, CloudPath
from cloudpathlib.exceptions import InvalidPrefixError


@register_client_class('hail-az')
class HailAzureBlobClient(AzureBlobClient):
    """
    Client is identical to the base Azure implementation.
    """

    pass


@register_path_class('hail-az')
class HailAzureBlobPath(AzureBlobPath):
    """
    Extending Path implementation to support hail-az:// and https:// schemes.
    >>> CloudPath('hail-az://myaccount/mycontainer/tmp')
    HailAzureBlobPath('hail-az://myaccount/mycontainer/tmp')
    >>> CloudPath('https://myaccount.blob.core.windows.net/mycontainer/tmp')
    HailAzureBlobPath('hail-az://myaccount/mycontainer/tmp')
    """

    cloud_prefix: str = 'hail-az://'
    client: 'HailAzureBlobClient'

    def __init__(
        self,
        cloud_path: str | CloudPath,
        client: AzureBlobClient | None = None,
        token: str | None = None,
    ):
        if isinstance(cloud_path, str):
            parsed = urlparse(cloud_path)
            m = re.match(
                r'(?P<account>[a-z0-9]+)(\.(?P<type>blob|dfs)(\.core\.windows\.net)?)?',
                parsed.netloc,
                flags=re.IGNORECASE,
            )
            if m is None:
                raise ValueError(f"Bad Azure path '{cloud_path}'")
            account = m.group('account')
            fstype = m.group('type') or 'blob'
            account_url = f'https://{account}.{fstype}.core.windows.net/'
            optional_type = '' if fstype == 'blob' else '.' + fstype
            cloud_path = f"{HailAzureBlobPath.cloud_prefix}{account}{optional_type}/{parsed.path.lstrip('/')}"
            if (
                client is None
                or parsed.query
                or token
                or client.service_client.account_name != account
            ):
                if token is not None:
                    token = '?' + token.lstrip('?')
                elif parsed.query:
                    token = '?' + parsed.query
                client = HailAzureBlobClient(account_url, token)

        super().__init__(cloud_path, client=client)

    @classmethod
    def is_valid_cloudpath(cls, path: str | CloudPath, raise_on_error=False) -> bool:
        """
        Also allowing HTTP.
        """
        valid = bool(
            re.match(
                fr'({HailAzureBlobPath.cloud_prefix}|https://[a-z0-9]+\.(blob|dfs)\.core\.windows\.net)',
                str(path).lower(),
            )
        )

        if raise_on_error and not valid:
            raise InvalidPrefixError(
                f'{path} is not a valid path since it does not start with {cls.cloud_prefix} '
                f'or valid Azure https blob or dfs location.'
            )

        return valid

    @property
    def container(self) -> str:
        """
        Minus the account part.
        """
        return self._no_prefix.split('/', 2)[1]

    @property
    def blob(self) -> str:
        """
        No prefix, no account part.
        """
        return super().blob.split('/', 1)[1]
