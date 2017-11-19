from os.path import exists

from .support.support.subprocess_ import run_command


def plot_fastq_gz_or_bam(fastq_gz_or_bam_file_path, overwrite=False):
    """
    Plot .fastq.gz or .bam file.
    Arguments:
        fastq_gz_or_bam_file_path (str):
        overwrite (bool):
    Returns:
        None
    """

    plot_zip_prefix_path = fastq_gz_or_bam_file_path + '.plot'
    plot_tsv_file_path = plot_zip_prefix_path + '.tsv'

    if not overwrite:
        if exists(plot_zip_prefix_path + '.zip'):
            raise FileExistsError(plot_zip_prefix_path + '.zip')
        if exists(plot_tsv_file_path):
            raise FileExistsError(plot_tsv_file_path)

    run_command(
        'fastqp --output {} --text {} {}'.format(plot_zip_prefix_path,
                                                 plot_tsv_file_path,
                                                 fastq_gz_or_bam_file_path),
        print_command=True)
