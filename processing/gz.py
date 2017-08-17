import gzip
from os import remove, rename

from Bio import bgzf

from .helper.helper.subprocess_ import run_command


def bgzip(file_path):
    """
    bgzip file_path.
    Arguments:
        file_path (str):
    Returns:
        None
    """

    command = 'bgzip -f {}'.format(file_path)

    run_command(command)


def tabix(bgzipped_file_path):
    """
    tabix bgzipped_file_path.
    Arguments:
        bgzipped_file_path (str):
    Returns:
        None
    """

    command = 'tabix -f {}'.format(bgzipped_file_path)

    run_command(command)


def convert_gzipped_to_bgzipped(gzipped_file_path):
    """
    Convert gzipped file to bgzipped file.
    Arguments:
        gzipped_file_path (str):
    Returns:
        None
    """

    tmp_file_path = '{}.convert_gzipped_to_bgzipped.tmp'.format(
        gzipped_file_path)

    # Read from gzipped file path
    with gzip.open(gzipped_file_path, 'rt') as gzipped_file:

        # Write to tmp file path
        with bgzf.open(tmp_file_path, 'wt') as tmp_bgzipped_file:
            for line in gzipped_file:
                tmp_bgzipped_file.write(line)

    # Remove gzipped file path
    remove(gzipped_file_path)

    # Rename tmp file path to gzipped file path
    rename(tmp_file_path, gzipped_file_path)
