from .support.support.subprocess_ import run_command


def bgzip(file_path):
    """
    bgzip file_path.
    Arguments:
        file_path (str):
    Returns:
        str:
    """

    command = 'bgzip -f {}'.format(file_path)

    run_command(command)

    return file_path + '.gz'


def tabix(bgzip_file_path):
    """
    tabix bgzip_file_path.
    Arguments:
        bgzip_file_path (str):
    Returns:
        str:
    """

    command = 'tabix -f {}'.format(bgzip_file_path)

    run_command(command)

    return bgzip_file_path + '.tbi'
