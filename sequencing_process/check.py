from inspect import stack
from os.path import dirname, exists, join

from .support.support.subprocess_ import run_command


def check_fastq_gz_using_fastqc(fastq_gz_file_path, n_jobs=1, overwrite=False):
    """
    Check .fastq.gz file using fastqp.
    Arguments:
        fastq_bam_file_path (str):
        n_jobs (int):
        overwrite (bool):
    Returns:
        None
    """

    run_command(
        'fastqc --threads {} {}'.format(n_jobs, fastq_gz_file_path),
        print_command=True)


def check_fastq_gz_or_bam_using_fastqp(fastq_gz_or_bam_file_path,
                                       overwrite=False):
    """
    Check .fastq.gz or .bam file using fastqp.
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
        'fastqp --output {} --text {} --count-duplicates {}'.format(
            plot_zip_prefix_path, plot_tsv_file_path,
            fastq_gz_or_bam_file_path),
        print_command=True)


def check_bam_using_samtools_flagstat(bam_file_path,
                                      n_jobs=1,
                                      output_file_path=None,
                                      overwrite=False):
    """
    Check .bam file using samtools flagstat.
    Arguments:
        bam_file_path (str):
        n_jobs (int):
        output_file_path (str):
        overwrite (bool):
    Returns:
        None
    """

    if not output_file_path:
        output_file_path = join(dirname(bam_file_path), stack()[0][3])

    if not overwrite and exists(output_file_path):
        raise FileExistsError(output_file_path)

    run_command(
        'samtools flagstat --threads {} {} > {}'.format(
            n_jobs, bam_file_path, output_file_path),
        print_command=True)
