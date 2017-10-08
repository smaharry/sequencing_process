from os.path import exists

from .process_bam import index_bam_using_samtools
from .support.support.subprocess_ import run_command


def align_fastq_gzs_using_hisat2(hisat2_index_file_path_prefix,
                                 fastq_gz_file_paths,
                                 sequence_type,
                                 n_jobs=1):
    """
    Align .fastq.gz files using hisat2.
    Arguments:
        hisat2_index_file_path_prefix (str): reference .fasta file path
        fastq_gz_file_paths (iterable): (< 2)
        sequence_type (str): 'DNA' | 'RNA'
        n_jobs (int):
    Returns:
        str:
    """

    if not all([
            exists('{}.{}.ht2'.format(hisat2_index_file_path_prefix, i))
            for i in range(1, 9)
    ]):
        print('Could not find HISAT2 index; building it ...')

        command = 'hisat2-build {} {}'.format(hisat2_index_file_path_prefix,
                                              hisat2_index_file_path_prefix)
        run_command(command)

    if len(fastq_gz_file_paths) == 1:
        print('Using single-end ...')
        sample_command = '-U {}'.format(*fastq_gz_file_paths)

    elif len(fastq_gz_file_paths) == 2:
        print('Using paired-end ...')
        sample_command = '-1 {} -2 {}'.format(*fastq_gz_file_paths)

    else:
        raise ValueError(
            'fastq_gz_file_paths must be an iterable containing unpaired or paired .fastq.gz file paths.'
        )

    if sequence_type == 'DNA':
        additional_arguments = '--dta-cufflinks'

    elif sequence_type == 'RNA':
        additional_arguments = '--no-spliced-alignment'

    output_bam_file_path = fastq_gz_file_paths[0] + '.align_fastq_gzs_using_hisat2.bam'

    command = 'hisat2 {} -x {} -p {} {} | samtools sort -@ {} > {}'.format(
        sample_command, hisat2_index_file_path_prefix, n_jobs,
        additional_arguments, n_jobs, output_bam_file_path)
    run_command(command)

    command = 'samtools index -@ {} {}'.format(n_jobs, output_bam_file_path)
    run_command(command)

    return output_bam_file_path


def align_fastq_gzs_using_bwa(fasta_file_path, fastq_gz_file_paths, n_jobs=1):
    """
    Align unpaired or paired .fastq.gz files using bwa.
    Arguments:
        fasta_file_path (str):
        fastq_gz_file_paths (iterable): (< 2)
        n_jobs (int):
    Returns:
        str:
    """

    if not all([
            exists('{}.{}'.format(fasta_file_path, suffix))
            for suffix in ['bwt', 'pac', 'ann', 'amb', 'sa']
    ]):
        print('Could not find BWA index; building it ...')
        command = 'bwa index {}'.format(fasta_file_path)
        run_command(command)

    output_bam_file_path = fastq_gz_file_paths[0] + '.align_fastq_gzs_using_bwa.bam'

    command = 'bwa mem -t {} {} {} | samtools sort -@ {} > {}'.format(
        n_jobs, fasta_file_path, ' '.join(fastq_gz_file_paths), n_jobs,
        output_bam_file_path)
    run_command(command)

    return index_bam_using_samtools(output_bam_file_path, n_jobs=n_jobs)
