import os
from functools import cache
from typing import Optional

from fsspec import AbstractFileSystem
from fsspec.implementations.local import LocalFileSystem
from s3fs import S3FileSystem


@cache
def get_bucket_filesystem(
    endpoint: Optional[str] = None,
    access_key: Optional[str] = None,
    secret_key: Optional[str] = None,
) -> AbstractFileSystem:
    if endpoint is None:
        endpoint = os.getenv('BUCKET_ENDPOINT', 'http://minio:9000')

    if access_key is None:
        access_key = os.getenv('BUCKET_ACCESS_KEY', None)

    if secret_key is None:
        secret_key = os.getenv('BUCKET_SECRET_KEY', None)

    return S3FileSystem(endpoint_url=endpoint, key=access_key, secret=secret_key)


@cache
def get_local_filesystem() -> AbstractFileSystem:
    return LocalFileSystem()
