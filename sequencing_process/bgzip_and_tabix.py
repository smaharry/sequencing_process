from .support.support.subprocess_ import run_command


def bgzip_and_tabix(file_path, n_job=1, overwrite=False):
    """
    Bgzip and tabix file_path.
    Arguments:
        file_path (str):
        n_job (int):
        overwrite (bool):
    Returns:
        str:
    """

    output_file_path = file_path + '.gz'

    bgzip_additional_arguments = []
    tabix_additional_arguments = []

    if overwrite:
        bgzip_additional_arguments.append('--force')
        tabix_additional_arguments.append('--force')

    run_command(
        'bgzip --threads {} {} {} && tabix {} {}'.format(
            n_job, ' '.join(bgzip_additional_arguments), file_path,
            ' '.join(tabix_additional_arguments), output_file_path),
        print_command=True)

    return output_file_path
