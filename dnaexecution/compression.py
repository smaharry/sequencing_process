import gzip
from os import remove, rename
from os.path import isfile

from Bio import bgzf

from .helper.helper.subprocess_ import run_command


def bgzip_tabix(file_path):
    """
    bgzip and tabix <file_path>.
    :param file_path: str
    :return: str
    """

    if file_path.endswith('.gz'):

        if not isfile('{}.tbi'.format(file_path)):  # Not tabixed
            tabix(file_path)

        return file_path

    else:

        bgzip_output_file_path = bgzip(file_path)
        tabix(bgzip_output_file_path)

    return bgzip_output_file_path


def bgzip(file_path):
    """
    bgzip <file_path>.
    :param file_path: str
    :return: str
    """

    command = 'bgzip -f {}'.format(file_path)

    run_command(command)

    return '{}.gz'.format(file_path)


def tabix(file_path):
    """
    tabix <file_path>.
    :param file_path: str
    :return: str
    """

    command = 'tabix -f {}'.format(file_path)

    run_command(command)

    return '{}.tbi'.format(file_path)


def convert_gzipped_to_bgzipped(gzipped_file_path):
    """
    Convert gzipped file to bgzipped file.
    :param gzipped_file_path: str
    :return: str
    """

    temp_file_path = '{}.convert_gzipped_to_bgzipped'.format(gzipped_file_path)

    with gzip.open(gzipped_file_path, 'rt') as gzipped_file:

        with bgzf.open(temp_file_path, 'wt') as temp_bgzipped_file:

            for line in gzipped_file:
                temp_bgzipped_file.write(line)

    remove(gzipped_file_path)

    rename(temp_file_path, gzipped_file_path)

    return gzipped_file_path
