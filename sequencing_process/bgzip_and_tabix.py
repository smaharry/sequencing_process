from .support.support.subprocess_ import run_command_and_monitor


def bgzip_and_tabix(file_path, n_jobs=1, overwrite=False):
    """
    Bgzip and tabix file_path.
    Arguments:
        file_path (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_file_path = file_path + '.gz'

    bgzip_additional_argument = ''
    tabix_additional_argument = ''
    if overwrite:
        bgzip_additional_argument += ' --force'
        tabix_additional_argument += ' --force'

    run_command_and_monitor(
        'bgzip --threads {} {} {}; tabix {} {}'.format(
            n_jobs, bgzip_additional_argument, file_path,
            tabix_additional_argument, output_file_path),
        print_command=True)

    return output_file_path
