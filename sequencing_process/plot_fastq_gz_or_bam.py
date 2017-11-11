from .support.support.subprocess_ import run_command_and_monitor


def plot_fastq_gz_or_bam(fastq_gz_or_bam_file_path, overwrite=False):
    """
    Plot .fastq.gz or .bam file.
    Arguments:
        fastq_gz_or_bam_file_path (str):
        overwrite (bool):
    Returns:
        None
    """

    run_command_and_monitor(
        'fastqp --output {0}.plot --text {0}.plot.tsv {0}'.format(
            fastq_gz_or_bam_file_path),
        print_command=True)
