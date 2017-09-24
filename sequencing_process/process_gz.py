from .support.support.subprocess_ import run_command


def bgzip(file_path):
    """
    bgzip file.
    Arguments:
        file_path (str):
    Returns:
        str:
    """

    command = 'bgzip -f {}'.format(file_path)

    run_command(command)

    return file_path + '.gz'


def tabix(bgzipped_file_path):
    """
    tabix bgzipped file.
    Arguments:
        bgzipped_file_path (str):
    Returns:
        str:
    """

    command = 'tabix -f {}'.format(bgzipped_file_path)

    run_command(command)

    return bgzipped_file_path + '.tbi'
