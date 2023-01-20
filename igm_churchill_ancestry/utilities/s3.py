import logging
from typing import Tuple
from urllib.parse import urlparse

import boto3
from awscli.clidriver import create_clidriver


def validate_s3_path(s3_url: str) -> Tuple[str, str]:
    """Validate an s3 path."""
    if not s3_url.startswith("s3://"):
        raise ValueError("Invalid s3 path.")
    parse = urlparse(s3_url)
    bucket, path = parse.netloc, parse.path.lstrip('/')
    return bucket, path


def download_s3_file(s3_url: str, local_file: str) -> None:
    """Download file from s3 locally."""
    bucket, key = validate_s3_path(s3_url)
    client = boto3.client("s3")
    logging.info(f"Downloading {s3_url} to {local_file}")
    client.download_file(bucket, key, local_file)


def upload_s3_file(local_file: str, s3_url: str) -> None:
    """Upload a local file to s3."""
    bucket, key = validate_s3_path(s3_url)
    client = boto3.client("s3")
    logging.info(f"Uploading {local_file} to {s3_url}")
    client.upload_file(local_file, bucket, key)


def download_s3_directory(s3_dir_name: str, local_dir_name: str) -> None:
    """Validate the s3 url, create cmd, download."""
    validate_s3_path(s3_dir_name)
    cmd = [
        's3', 'cp', '--recursive', '--no-progress',
        s3_dir_name, local_dir_name
    ]
    logging.info(f"Downloading {s3_dir_name} to {local_dir_name}")
    rc = create_clidriver().main(cmd)
    if rc != 0:
        raise RuntimeError(
            f"There was an error downloading the directory. {cmd}")


def upload_s3_directory(local_dir_name: str, s3_dir_name: str) -> None:
    """Validate the s3 url, create cmd, upload."""
    validate_s3_path(s3_dir_name)
    cmd = [
        's3', 'cp', '--recursive', '--no-progress',
        local_dir_name, s3_dir_name
    ]
    logging.info(f"Uploading {local_dir_name} to {s3_dir_name}")
    rc = create_clidriver().main(cmd)
    if rc != 0:
        raise RuntimeError(
            f"There was an error uploading the directory. {cmd}")