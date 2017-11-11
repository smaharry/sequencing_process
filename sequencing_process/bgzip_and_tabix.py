from .support.support.subprocess_ import run_command_and_monitor


def bgzip_and_tabix(file_path, n_jobs=1):
    """
    Bgzip and tabix file_path.
    Arguments:
        file_path (str):
        n_jobs (int):
    Returns:
        str:
    """

    output_file_path = file_path + '.gz'

    run_command_and_monitor(
        'bgzip --force --threads {} {}; tabix --force {}'.format(
            n_jobs, file_path, output_file_path),
        print_command=True)

    return output_file_path
