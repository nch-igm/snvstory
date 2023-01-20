import os
import shutil
import hashlib
import logging

from .s3 import (upload_s3_file, upload_s3_directory,
                 download_s3_file, download_s3_directory)


def verbose_return(old_path: str, new_path: str) -> str:
    logging.info(f'\t{old_path} -> {new_path}')
    return new_path


def _s3_parent_name_md5sum(in_path: str) -> str:
    basename = os.path.basename(in_path)
    parent = in_path.replace(basename, '')
    # https://stackoverflow.com/questions/5297448/
    return hashlib.md5(parent.encode('utf-8')).hexdigest()


def flex_input(
        in_path: str,
        out_dir: str = None,
        directory: bool = False,
        force_copy: bool = False,
        prepend_hash: bool = False,
        basename: str = None) -> str:
    """
    Allows specifying a path to either an S3 file or a local file. If valid
    local file, returns the same path as input. If S3 file, downloads the file
    and returns the new local path. If input is a directory, must set
    'directory' parameter to True.

    Args
    ----
    in_path: The path to an s3 file or a local file.
    out_dir: Alternate output path for downloaded or copied file.
                   If None, than the file is downloaded to the current
                   directory '.'
    directory: Whether in_path is a directory or not.
    force_copy: Copy in_path to out_dir using the same basename.
    prepend_hash: Prepend a hash of the S3 "parent" path to file.
                         This should prevent basename collisions, such as for
                         FASTQs from multiple sequencers for the same sample.
    basename: New name of file or directory in out_dir

    Returns
    -------
    out_path (str): The location of the downloaded s3 file or copied local
    file. If out_dir not supplied and in_path already exists, out_path will
    equal in_path.
    """
    has_s3_prefix = in_path.startswith("s3://")
    # validate params
    out_dir = out_dir if out_dir else "."
    if directory and not (os.path.isdir(in_path) or has_s3_prefix):
        raise ValueError(f"{in_path} is not a directory.")
    if force_copy and out_dir == os.path.dirname(in_path):
        raise ValueError(f"Cannot force_copy since in_path {in_path} is already in out_dir {out_dir}.")
    if not (has_s3_prefix or os.path.exists(in_path)):
        raise ValueError(f"in_path {in_path} is not s3 or a local file.")
    # create basename
    if not basename:
        basename = os.path.basename(in_path)
    if prepend_hash:
        s3_parent_hash = _s3_parent_name_md5sum(in_path)
        basename = f"{s3_parent_hash}_{basename}"
    # determine out_path
    if os.path.exists(in_path) and not force_copy:
        out_path = in_path
    else:
        out_path = os.path.join(out_dir, basename)
    # handle s3 paths
    if has_s3_prefix:
        if directory:
            download_s3_directory(in_path, out_path)
        else:
            download_s3_file(in_path, out_path)
    # handle local files
    elif os.path.exists(in_path) and force_copy:
        if directory:
            shutil.copytree(in_path, out_dir)
        else:
            shutil.copy2(in_path, out_path)
    return verbose_return(in_path, out_path)


def flex_output(in_path: str, out_dir: str = None,
                file_rename: str = None) -> str:
    """
    Copies a file or directory to either a local output directory or an S3
    path. If it's valid local path, returns the same path as input.
    If S3 path, uploads the input file or directory to that path.

    Args
    ----
    in_path (str): The path to a local file.
    out_dir (str): Local directory or path in S3.
    file_rename (str): The new basename to be used for the file.

    Returns
    -------
    out_path (str): The location of the uploaded s3 file or local file.

    Examples:

    ```python
    >>> flex_output("./file.txt", "s3://test/path/")
    ```
    """
    # validate params
    if not os.path.exists(in_path):
        raise ValueError(f"in_path does not exist: {in_path}")
    if not out_dir:
        out_dir = "."
    if not file_rename:
        file_rename = os.path.basename(in_path)
    out_path = os.path.join(out_dir, file_rename)
    if in_path == out_path:
        print("Source and destination are the same. No copy needed.")
        return verbose_return(in_path, out_path)
    # handle s3
    if out_dir.startswith("s3://"):
        if os.path.isdir(in_path):
            upload_s3_directory(in_path, out_path)
        elif os.path.isfile(in_path):
            upload_s3_file(in_path, out_path)
        else:
            raise ValueError(
                f"Input {in_path} is not a valid file or directory path.")
    # handle local file
    if os.path.isdir(out_dir):
        if os.path.isdir(in_path):
            shutil.copytree(in_path, out_path)
        else:
            shutil.copy(in_path, out_path)

    return verbose_return(in_path, out_path)


class FlexIoMixin:

    @staticmethod
    def flex_input(*args, **kwargs) -> str:
        return flex_input(*args, **kwargs)

    @staticmethod
    def flex_output(*args, **kwargs) -> str:
        return flex_output(*args, **kwargs)
