from os.path import isfile

from .support.support.subprocess_ import run_command_and_monitor


def plot_fastq_gzs_or_bams(fastq_gz_or_bam_file_paths, overwrite=False):
    """
    Plot .fastq.gz or .bam files.
    Arguments:
        fastq_gz_or_bam_file_paths (iterable):
        overwrite (bool):
    Returns:
        str:
    """

    output_path = fastq_gz_or_bam_file_paths[0] + '.plot'

    if isfile(output_path) and not overwrite:
        raise FileExistsError('{} exists.'.format(output_path))

    run_command_and_monitor(
        'fastqp --output {0} --text {0} {1}'.format(
            output_path, ' '.join(fastq_gz_or_bam_file_paths)),
        print_command=True)

    return output_path
