from os.path import isfile

from .support.support.subprocess_ import run_command_and_monitor


def plot_fastq_gz_or_bam(fastq_gz_or_bam_file_path, overwrite=False):
    """
    Plot .fastq.gz or .bam file.
    Arguments:
        fastq_gz_or_bam_file_path (str):
        overwrite (bool):
    Returns:
        str:
    """

    output_path = fastq_gz_or_bam_file_path + '.plot'

    if isfile(output_path) and not overwrite:
        raise FileExistsError('{} exists.'.format(output_path))

    run_command_and_monitor(
        'fastqp --output {0} --text {0} {1}'.format(output_path,
                                                    fastq_gz_or_bam_file_path),
        print_command=True)

    return output_path
