from os.path import isfile

from ._print_and_run_command import _print_and_run_command


def bgzip_and_tabix(
        file_path,
        n_job=1,
        overwrite=False,
):

    output_file_path = file_path + '.gz'

    if not overwrite and isfile(output_file_path):

        raise FileExistsError(output_file_path)

    if overwrite:

        force = '--force'

    else:

        force = ''

    _print_and_run_command(
        'bgzip --threads {0} {1} {2} && tabix {1} {3}'.format(
            n_job,
            force,
            file_path,
            output_file_path,
        ))

    return output_file_path
