from .support.support.subprocess_ import run_command


def bgzip_and_tabix(file_path):
    """
    Bgzip and tabix file_path.
    Arguments:
        file_path (str):
    Returns:
        str:
    """

    output_file_path = file_path + '.gz'

    run_command('bgzip -f {}; tabix -f {}'.format(file_path, output_file_path))

    return output_file_path
